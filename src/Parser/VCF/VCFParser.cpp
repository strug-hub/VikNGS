#include "VCFParserUtils.h"
#include <future>
#include <thread>
#include <chrono>

static const std::string VCF_PARSER = "VCF parser";

std::vector<Variant> runBatch(TestInput &input, Request &req, std::vector<Variant> &variants){

    std::vector<std::vector<int>> collapse;

    if(req.shouldCollapseK())
        collapse = collapseEveryK(req.collapse, variants.size());
    else if(req.shouldCollapseBed())
        collapse = collapseEveryK(req.collapse, variants.size());

    TestInput in = addVariants(input, variants, collapse);
    std::vector<Variant> results = runTest(in, req);

    return results;
}

class ParallelProcess
{
    int index;
    bool started;
    std::future<std::vector<Variant>> futureResults;

    TestInput input;
    Request req;
    std::vector<Variant> variants;

public:

    ParallelProcess(bool s) : started(s) { }

    inline bool isDone(){
        return !started || futureResults.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }

    inline bool hasStarted(){
        return started;
    }

    inline std::vector<Variant> getResults(){
        return futureResults.get();
    }

    inline void reset(){
        started = false;
    }

    void copyValues(int index, TestInput input, Request req, std::vector<Variant> variants){
        this->index = index;
        this->input = input;
        this->req = req;
        this->variants = variants;
    }

    void begin(){
        started = true;
        this->index = index;

        futureResults = std::async(std::launch::async,
                [this] { return runBatch(input, req, variants); }) ;
    }
};

int findInterval(std::vector<Interval> &intervals, Variant &variant){

    //todo: optimize!@!!!!!
    for(int i = 0; i < intervals.size(); i++){

        if(intervals[i].start <= variant.pos && intervals[i].end >= variant.pos
                && intervals[i].chr == variant.chr)
            return i;
    }
    return -1;
}

std::vector<Variant> parseVCFLines(TestInput &input, Request &req){

    std::string vcfDir = req.vcfDir;
    int minPos = req.minPos;
    int maxPos = req.maxPos;
    int nthreads = req.nthreads - 1;

	File vcf;
	vcf.open(vcfDir);
    std::vector<Variant> variants;
    std::vector<Variant> output;
    std::vector<std::string> filterInfo;
    std::vector<int> filterCode;
    std::vector<ParallelProcess> threads;
    for(int i =0; i< nthreads; i++)
        threads.emplace_back(false);

	//skips header
	extractHeader(vcf);

	int lineCount = 0;

    //skip until POS in VCF >= minPos
	if(minPos > -1){
		while (vcf.hasNext()) {
			std::string line = vcf.nextLine();
            std::vector<std::string> columns = split(line, VCF_SEPARATOR);
		
        int pos = std::stoi(columns[LOC]);
        if(pos >= minPos)
            break;
		}
	}
	
	while (vcf.hasNext()) {

        if(lineCount % 5000 == 0){
            if(lineCount == 0)
                printInfo("Parsing VCF file...");
            else
                printInfo( std::to_string(lineCount) + " variant lines have been parsed so far");
        }

        std::string line = vcf.nextLine();

        std::vector<std::string> columns = split(line, VCF_SEPARATOR);
        variants.emplace_back(constructVariant(columns, req.regularTest));

        if(!variants.back().isValid()){
            std::string lineNumber = std::to_string(vcf.getLineNumber());
            printWarning(VCF_PARSER, "Skipping line " + lineNumber + " of VCF file - " + variants.back().getErrorMessage());
            variants.pop_back();
            continue;
        }

        if(maxPos > -1 && variants.back().pos > maxPos)
            break;

        //todo move down
        lineCount++;

        if(req.shouldCollapseBed()){
            int index = findInterval(input.intervals, variants.back());
            if(index < 0){
                variants.pop_back();
                continue;
            }
            else
                variants.back().interval = input.intervals[index];

        }

        calculateExpectedGenotypes(variants.back());
        int code = filterVariant(req, variants.back(), input.Y, input.family);

        if(code > 0){
            filterInfo.emplace_back(variants.back().toString());
            filterCode.emplace_back(code);
            variants.pop_back();
        }
        else{

            if(variants.size() >= req.batchSize){

                for(int m = 0; m < nthreads; m++){
                    if(threads[m].hasStarted() && threads[m].isDone()){
                        std::vector<Variant> results = threads[m].getResults();
                        outputPvals(results, req.outputDir);
                        threads[m].reset();

                        if(req.retainVariants){
                            for(int l = 0; l<results.size(); l++){
                                results[l].reduceSize();
                                output.push_back(results[l]);
                            }
                        }
                    }
                }

                //TODO: possible race condition with above loop?
                //wait for available thread
                bool printWait =false;
                while(true){

                    bool stopLooping = false;
                    for(int m = 0; m < nthreads; m++){
                        if(threads[m].isDone()){

                            int total = filterInfo.size() + variants.size();

                            printInfo("====================================");
                            printInfo("Read variants " + std::to_string(lineCount-total + 1) + " to " + std::to_string(lineCount));
                            printFilterResults(req, filterInfo, filterCode, total);
                            printInfo("Starting batch of tests on thread " + std::to_string(m));
                            printInfo("====================================");

                            threads[m].reset();
                            threads[m].copyValues(m, input, req, variants);
                            threads[m].begin();

                            variants.clear();
                            filterInfo.clear();
                            filterCode.clear();

                            stopLooping = true;
                            break;
                        }
                    }

                    if(stopLooping)
                        break;

                    if(!printWait)
                        printInfo("Waiting for available thread to start next batch of tests");

                    printWait = true;
                    std::this_thread::sleep_for(std::chrono::milliseconds(500));
                }
            }

        }
	}

    printInfo( std::to_string(lineCount) + " variants read from VCF");
	vcf.close();

    int total = filterInfo.size() + variants.size();
    printInfo("====================================");
    printInfo("Read variants " + std::to_string(lineCount-total + 1) + " to " + std::to_string(lineCount));
    printFilterResults(req, filterInfo, filterCode, total);
    printInfo("Starting batch of tests");
    printInfo("====================================");
    printInfo("VCF parsing complete");

    std::vector<Variant> results = runBatch(input, req, variants);

    variants.clear();
    filterInfo.clear();
    filterCode.clear();

    if(req.retainVariants){
        for(int l = 0; l< results.size(); l++){
            results[l].reduceSize();
            output.push_back(results[l]);
        }
    }

    //catch all remaining jobs
    bool printWait = false;
    while(true){

        bool stopLooping = true;
        for(int m = 0; m < nthreads; m++){
            if(!threads[m].isDone())
                stopLooping = false;
            else if (threads[m].hasStarted() && threads[m].isDone()){

                std::vector<Variant> results = threads[m].getResults();
                outputPvals(results, req.outputDir);
                threads[m].reset();

                if(req.retainVariants){
                    for(int l = 0; l< results.size(); l++){
                        results[l].reduceSize();
                        output.push_back(results[l]);
                    }
                }
            }
        }

        if(stopLooping)
            break;

        if(!printWait)
            printInfo("Waiting for remaining batches to finish");

        printWait = true;
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }

    return output;
}


