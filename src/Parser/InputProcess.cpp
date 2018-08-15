#include "Parser.h"
#include "Filter.h"

#include <future>
#include <thread>
#include <chrono>
#include <queue>

static const std::string ERROR_SOURCE = "INPUT_PARSER";

SampleInfo parseSampleInfo(Request & req) {
	
    SampleInfo sampleInfo;
    std::string dir = req.getSampleDir();

    std::map<std::string, int> IDmap = getSampleIDMap(req.getVCFDir());
    validateSampleIDs(dir, IDmap);

    sampleInfo.setY(parseSamplePhenotype(dir, IDmap));
    sampleInfo.setG(parseSampleGroupID(dir, IDmap));

    VectorXd G = sampleInfo.getG();
    sampleInfo.setGroupDepthMap(
                parseSampleReadDepth(dir, IDmap, G, req.getHighLowCutOff() )
                );
    sampleInfo.setZ(parseSampleCovariates(dir, IDmap));

    //todo: print more sample info?

    if(sampleInfo.hasCovariates())
        printInfo(std::to_string(sampleInfo.ncov()) + " covariates parsed");

    return sampleInfo;
}

std::vector<Variant> constructVariants(SampleInfo * sampleInfo, Request * req, std::vector<std::string> &lines){

    bool getVCFCalls = req->requireVCFCalls();
    bool calculateExpected = req->requireExpectedGenotypes();
    bool calculateCalls = req->requireGenotypeCalls();

    std::vector<Variant> variants;
    variants.reserve(lines.size());

    //contruct variants
    for(size_t i = 0; i < lines.size(); i++){

        //extract to the FILTER column
        std::vector<std::string> info = splitString(lines[i], VCF_SEP, FILTER);
        Filter filter = filterByVariantInfo(req, info[CHROM], info[POS], info[REF], info[ALT], info[FILTER]);

        if(filter == Filter::IGNORE || (filter != Filter::VALID && !req->shouldKeepFiltered()))
            continue;

        std::vector<std::string> columns = splitString(lines[i], VCF_SEP);
        Variant variant = constructVariant(columns, calculateExpected, calculateCalls, getVCFCalls);
        variant.setFilter(filter);

        if(variant.isValid()){

            VectorXd Y = sampleInfo->getY();
            filter = filterByGenotypes(req, variant, Y, sampleInfo->getFamily());

            if(filter != Filter::VALID && !req->shouldKeepFiltered())
                continue;

            variant.setFilter(filter);
        }
        variants.push_back(variant);
    }
    variants.shrink_to_fit();
    return variants;
}

int findInterval(IntervalSet * is, std::string chr, int pos, int searchHint){
    size_t hint = static_cast<size_t>(searchHint);
    std::vector<Interval>* intervals = is->get(chr);
    size_t left = 0; size_t right = intervals->size();
    size_t middle = hint;
    if(hint <= left || hint >= right)
        middle = (left + right) / 2;

    while (left <= right) {
        if (intervals->at(middle).isIn(pos))
              return static_cast<int>(middle);
        else if (intervals->at(middle).isSmaller(pos))
              right = middle - 1;
        else
              left = middle + 1;

        middle = (left + right) / 2;
    }

    return -1;
}

std::vector<VariantSet> collapseVariants(Request * req, std::vector<Variant> variants, VariantSet leftover){

    IntervalSet * is = req->getIntervals();
    CollapseType collapse = req->getCollapseType();
    std::vector<VariantSet> variantSet;

    if(collapse == CollapseType::NONE){
        if(leftover.size() > 0)
            variantSet.push_back(leftover);
        for(size_t i = 0; i < variants.size(); i++)
            variantSet.emplace_back(VariantSet(variants[i]));
    }

    else if(collapse == CollapseType::COLLAPSE_K){
        int k = req->getCollapseSize();

        size_t startFrom = 0;
        while(leftover.validSize() < k && startFrom < variants.size()){
            leftover.addVariant(variants[startFrom]);
            startFrom++;
        }
        variantSet.push_back(leftover);

        int counter = 0;
        for(size_t i = startFrom; i < variants.size(); i++){
            if(counter == 0)
                variantSet.emplace_back(VariantSet(variants[i]));
            else if(counter < k)
                variantSet.back().addVariant(variants[i]);

            counter++;
            if(counter == k) counter = 0;
        }
    }
    else if(collapse == CollapseType::COLLAPSE_EXON || collapse == CollapseType::COLLAPSE_GENE){

        size_t startFrom = 0;
        while(startFrom < variants.size()){
            if(leftover.isIn(variants[startFrom]))
                leftover.addVariant(variants[startFrom]);
            else
                break;
            startFrom++;
        }
        if(leftover.size() > 0)
            variantSet.push_back(leftover);

        int hint = -1;
        for(size_t i = startFrom; i < variants.size(); i++){

            if(variantSet.back().isIn(variants[i]))
                variantSet.back().addVariant(variants[i]);
            else{
                int index = findInterval(req->getIntervals(), variants[i].getChromosome(), variants[i].getPosition(), hint);
                VariantSet newSet(variants[i]);
                newSet.setInterval(&req->getIntervals()->get(variants[i].getChromosome())->at(index));
                hint = index+1;
            }
        }
    }

    return variantSet;
}


std::vector<Variant> testBatch(SampleInfo input, Request req, std::vector<VariantSet> &variants){
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
private:
    Request* req;
    SampleInfo* sampleInfo;

    std::vector<std::string> lines;
    std::vector<Variant> variants;
    VariantSet leftover;

    bool running;
    std::future<std::vector<Variant>> futureVariants;
    std::future<std::vector<VariantSet>> futureVariantSets;

public:

    ParallelProcess(Request *r, SampleInfo *si) : req(r), sampleInfo(si) { running = false;}

    inline bool isRunning(){
        return running;
    }

    void parseAndFilter(std::vector<std::string> & lns){
        running = true;
        this->lines = lns;

        futureVariants = std::async(std::launch::async,
                [this] { return constructVariants(sampleInfo, req, lines); }) ;
    }
    inline bool isParseDone(){
        return running && futureVariants.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }

    inline std::vector<Variant> getFilterResults(){
        running = false;
        return futureVariants.get();
        lines.clear();
    }

    void createVariantSet(std::vector<Variant> & v, VariantSet & leftovers){
        running = true;
        this->variants = v;
        this->leftover = leftovers;
        futureVariantSets = std::async(std::launch::async,
                [this] { return collapseVariants(req, variants, leftover); }) ;
    }
    inline bool isCollapseDone(){
        return running && futureVariantSets.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }
    inline std::vector<VariantSet> getCollapseResults(){
        running = false;
        return futureVariantSets.get();
        variants.clear();
    }

    void beginAssociationTest(std::vector<VariantSet> variants){
        running = true;
        this->variants = variants;

        futureVariantSets = std::async(std::launch::async,
                [this] { return testBatch(sampleInfo, req, lines, startVariantNumber, startLineNumber); }) ;
    }
};

std::vector<VariantSet> processVCF(Request &req, SampleInfo &input) {

    size_t nthreads = static_cast<size_t>(req.getNumberThreads() - 1);

    File vcf;
    vcf.open(req.getVCFDir());

    size_t totalLineCount = 0;
    size_t batchSize = static_cast<size_t>(req.getBatchSize());
    std::vector<std::string> lines;

    std::vector<ParallelProcess> threads;
    for(size_t i = 0; i < nthreads; i++)
        threads.emplace_back(input, req);
    std::queue<ParallelProcess*> filterOrder;

    std::vector<Variant> constructedVariants;
    std::vector<VariantSet> collapsedVariants;
    VariantSet leftover;
    bool readyToCollapse = false;

    //skips header
    extractHeaderLine(vcf);

    while (vcf.hasNext() || lines.size() > 0 || constructedVariants.size() > 0){

        //read in line if batch not full
        if(lines.size() < batchSize){
            lines.emplace_back(vcf.nextLine());
            totalLineCount++;

            if(totalLineCount % 10000 == 0){
                if(totalLineCount == 0)
                    printInfo("Parsing VCF file...");
                else
                    printInfo(std::to_string(totalLineCount) + " variant lines have been parsed so far");
            }
        }

        //check if threads are done
        if(filterOrder.size() > 0 && filterOrder.front()->isParseDone()){
            std::vector<Variant> v = filterOrder.front()->getFilterResults();
            filterOrder.pop();
            constructedVariants.insert(constructedVariants.end(), v.begin(), v.end());
        }

        for(size_t m = 0; m < nthreads; m++){
            if(threads[m].isCollapseDone()){
                std::vector<VariantSet> v = threads[m].getCollapseResults();
                collapsedVariants.insert(collapsedVariants.end(), v.begin(), v.end());

                leftover = collapsedVariants.back();
                collapsedVariants.pop_back();
                readyToCollapse = true;
            }
        }

        //finished reading in a batch
        if(lines.size() >= batchSize || !vcf.hasNext()){
            //start new thread
            for(size_t m = 0; m < nthreads; m++){
                 if(!threads[m].isRunning()){
                     threads[m].parseAndFilter(lines);
                     filterOrder.push(&threads[m]);
                     lines.clear();
                     break;
                 }
             }
        }

        if(readyToCollapse && (constructedVariants.size() >= batchSize || !vcf.hasNext())){
            for(size_t m = 0; m < nthreads; m++){
                 if(!threads[m].isRunning()){
                     threads[m].createVariantSet(constructedVariants, leftover);
                     filterOrder.push(&threads[m]);
                     constructedVariants.clear();
                     break;
                 }
             }
        }


    printInfo( std::to_string(totalLineCount) + " variants identified in VCF file");
    vcf.close();

    bool printWait = true;
    while(true){

        bool stopLooping = true;
        int threadCount = 0;
        for(size_t m = 0; m < nthreads; m++){
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

