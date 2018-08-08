#include "InputParser.h"
#include <future>
#include <thread>
#include <chrono>

static const std::string ERROR_SOURCE = "INPUT_PARSER";

/*
Parses all input files and filters variants based on given criteria

Parameters contained in Request
------------------------
Parsing params:
@param vcfDir Directory of VCF file.
@param sampleDir Directory of sample directory file.
@param bedDir Directory of BED file.
@param highLowCutOff Cut-off for read depth. Anything >= value is considered high read depth.

Filtering params:
@param missingThreshold Proportion of sample data missing for variant to be filtered.
@param onlySNPs If true, indels will be filtered based on REF and ALT columns.
@param mustPASS If true, variants where FILTER is not 'PASS' will be filtered.
@param mafCut The minor allele frequency cut-off for common or rare variants.
@param common Indicates whether to keep common or rare variants, (common = true, rare = false).
------------------------

Output params:
@param X Matrix of explanatory variable.
@param Y Vector of response variable.
@param Z Matrix of covariates.
@param G Vector of group ID.
@param P Vector with probability of 0, 1 or 2 minor alleles.
@param interval The indexes of variants within each interval.

@return TestInput object holding all the output parameters.
*/
SampleInfo parseSampleInfo(Request req) {
	
	VectorXd Y, G; MatrixXd Z;
	std::map<int, int> readGroup;

    std::map<std::string, int> IDmap = getSampleIDMap(req.vcfDir);
	parseSampleLines(req, IDmap, Y, Z, G, readGroup);

    std::string family = determineFamily(Y);
	
    std::vector<Interval> intervals;
	if (req.shouldCollapseBed())
        intervals = parseBEDLines(req.bedDir, req.shouldCollapseExon());

    return buildTestInput(Y, Z, G, readGroup, intervals, family);
}

std::vector<Variant> runBatch(SampleInfo input, Request req, std::vector<std::string> &lines, int startVariantNumber, int startLineNumber){

    std::vector<Variant> variants;
    std::vector<std::string> filterInfo;
    std::vector<int> filterCode;
    int total = lines.size();

    //contruct variants
    for(int i = 0; i < total; i++){

        std::vector<std::string> columns = split(lines[i], VCF_SEPARATOR);

        variants.emplace_back(constructVariant(columns));
        if(!variants.back().isValid()){
            std::string lineNumber = std::to_string(startLineNumber + i);
            printWarning(INPUT_PARSER, "Skipping line " + lineNumber + " of VCF file - " + variants.back().getErrorMessage());
            variants.pop_back();
            continue;
        }

        calculateExpectedGenotypes(variants.back());
        variants.back().genotypeCalls = calculateGTCalls(variants.back().likelihood, variants.back().P);

        int code = filterVariant(req, variants.back(), input.Y, input.family);
        if(code > 0){
            filterInfo.emplace_back(variants.back().toString());
            filterCode.emplace_back(code);
            variants.pop_back();
        }
  //      else
  //          outputDebug(lines[i],"/home/scott/vikNGS/build-GUI-Desktop_Qt_5_11_0_GCC_64bit-Release");
    }

    lines.clear();

    printInfo("====================================");
    printInfo("Read variants " + std::to_string(startVariantNumber) + " to " + std::to_string(startVariantNumber + total));
    printFilterResults(req, filterInfo, filterCode, total);
    printInfo("Starting batch of tests");
    printInfo("====================================");

    if(variants.size() <= 0)
        return variants;

    std::vector<std::vector<int>> collapse;

    if(req.shouldCollapseK())
        collapse = collapseEveryK(req.collapse, variants.size());
    else if(req.shouldCollapseBed())
        collapse = collapseEveryK(req.collapse, variants.size());
    //todo

    SampleInfo in = addVariants(input, variants, collapse);

    //todo: remove
    //printInfo("rvs = " + std::to_string(req.rvs) + " --- useGT = " + std::to_string(req.regularTest));

    std::vector<Variant> results = runTest(in, req);
    outputPvals(results, req.outputDir);

    return results;
}

class ParallelProcess
{
    SampleInfo input;
    Request req;

    int startVariantNumber;
    int startLineNumber;
    std::vector<std::string> lines;

    bool initialized;
    std::future<std::vector<Variant>> futureResults;

public:

    ParallelProcess(SampleInfo ti, Request r) : input(ti), req(r) { initialized = false;}

    inline bool isDone(){
        return initialized && futureResults.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }

    inline bool isInitialized(){
        return initialized;
    }

    inline std::vector<Variant> getResults(){
        initialized = false;
        return futureResults.get();
    }

    void copyValues(SampleInfo input, Request req){
        this->input = input;
        this->req = req;
    }

    void begin(int startVariantNum, int startLineNum, std::vector<std::string> lns){
        initialized = true;
        this->lines = lns;
        this->startLineNumber = startLineNum;
        this->startVariantNumber = startVariantNum;

        futureResults = std::async(std::launch::async,
                [this] { return runBatch(input, req, lines, startVariantNumber, startLineNumber); }) ;
    }
};

std::vector<Variant> processVCF(SampleInfo input, Request req) {

    std::string vcfDir = req.vcfDir;
    std::string filterChr = req.filterChr;
    int minPos = req.minPos;
    int maxPos = req.maxPos;
    int nthreads = req.nthreads - 1;
    int waitTime = 500; //miliseconds

    File vcf;
    vcf.open(vcfDir);

    int lineCount = 0;
    int firstLineInBatch = 0;
    int batchCout = 0;
    std::vector<std::string> lines;
    std::vector<Variant> output;

    std::vector<ParallelProcess> threads;
    for(int i =0; i< nthreads; i++)
        threads.emplace_back(input, req);

    //skips header
    extractHeaderLine(vcf);


   // int maxVariants = 20000;
 //   int counter= 0;
    while (vcf.hasNext() ){ // && lineCount < maxVariants) {

   //     counter++;

  //     if(counter > maxVariants)
  //          break;

        if(lineCount % 10000 == 0){
            if(lineCount == 0)
                printInfo("Parsing VCF file...");
            else
                printInfo(std::to_string(lineCount) + " variant lines have been parsed so far");
        }

        lines.emplace_back(vcf.nextLine());
        lineCount++;

        try{
            if(!isIn(lines.back(), minPos, maxPos, filterChr)){
                lines.pop_back();
                continue;
            }

        } catch(...){
                printWarning(INPUT_PARSER, "Error while parsing line " + std::to_string(vcf.lineNumber) +
                             " of VCF file. Ensure CHR and POS columns are formatted correctly. Skipping variant.");
        }

        if(req.shouldCollapseBed()){

            /*int index = findInterval(input.intervals, variants.back());
            if(index < 0){
                variants.pop_back();
                continue;
            }
            else
                variants.back().interval = input.intervals[index];
*///todo
        }

        if(lines.size() >= req.batchSize || !vcf.hasNext()){

            bool printWait = true;
            while(true){

                bool stopLooping = false;

                //check if threads are done
                for(int m = 0; m < nthreads; m++){
                    if(threads[m].isDone()){
                        std::vector<Variant> results = threads[m].getResults();

                        if(req.retainVariants)
                            for(int l = 0; l < results.size(); l++){
                                output.push_back(results[l]);

                        }
                    }
                }

                //start new thread
               for(int m = 0; m < nthreads; m++){
                    if(!threads[m].isInitialized()){
                        threads[m].begin(firstLineInBatch, vcf.lineNumber-lines.size(), lines);
                        batchCout++;
                        firstLineInBatch = lineCount+1;
                        lines.clear();
                        stopLooping = true;
                        break;
                    }
                }

                if(stopLooping)
                    break;

                if(printWait)
                    printInfo("Waiting for available thread to start next batch of tests");

                printWait = false;
                std::this_thread::sleep_for(std::chrono::milliseconds(waitTime));
            }
        }
    }

    printInfo( std::to_string(lineCount) + " variants identified in VCF file");
    vcf.close();

    bool printWait = true;
    while(true){

        bool stopLooping = true;
        int threadCount = 0;
        for(int m = 0; m < nthreads; m++){
            if(threads[m].isDone()){
                std::vector<Variant> results = threads[m].getResults();

                if(req.retainVariants){
                    for(int l = 0; l < results.size(); l++)
                        output.push_back(results[l]);

                }
                printWait=true;
            }
            else if(threads[m].isInitialized()){
                stopLooping = false;
                threadCount++;
            }
        }

        if(stopLooping)
            break;

        if(printWait)
            printInfo("Waiting for " + std::to_string(threadCount) + " threads to finish");

        printWait = false;
        std::this_thread::sleep_for(std::chrono::milliseconds(waitTime));

    }
    return output;
}

