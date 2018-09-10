#include "Parser.h"
#include "Filter.h"
#include "../Test/Test.h"

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

    VectorXi G = sampleInfo.getG();
    sampleInfo.setGroupDepthMap(
                parseSampleReadDepth(dir, IDmap, G, req.getHighLowCutOff() )
                );
    sampleInfo.setZ(parseSampleCovariates(dir, IDmap));

    //todo: print more sample info?

    if(sampleInfo.hasCovariates())
        printInfo(std::to_string(sampleInfo.ncov()) + " covariates parsed");

    return sampleInfo;
}

std::vector<Variant> constructVariants(Request* req, SampleInfo* sampleInfo, std::vector<std::string> &lines){

    bool getVCFCalls = req->requireVCFCalls();
    bool calculateExpected = req->requireExpectedGenotypes();
    bool calculateCalls = req->requireGenotypeCalls();

    std::vector<Variant> variants;
    variants.reserve(lines.size());

    //contruct variants
    for(size_t i = 0; i < lines.size(); i++){

        //extract to the FILTER column
        std::vector<std::string> info = splitString(lines[i], VCF_SEP, FILTER + 1);
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

std::vector<VariantSet> collapseVariants(Request* req, std::vector<Variant>& variants, VariantSet& leftover){

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
                int index = findInterval(is, variants[i].getChromosome(), variants[i].getPosition(), hint);
                VariantSet newSet(variants[i]);
                newSet.setInterval(&is->get(variants[i].getChromosome())->at(index));
                hint = index+1;
            }
        }
    }

    return variantSet;
}


bool testBatch(Request* req, SampleInfo* sampleInfo, std::vector<VariantSet*> &variants){
    if(variants.size() <= 0)
        return true;

    int nboot = 0;
    if(req->useBootstrap())
        nboot = req->bootstrapSize();

    for(VariantSet* vs : variants)
        for(Test t : req->getTests()){
            if(vs->validSize() > 0){
                double pval = runTest(sampleInfo, vs, t, nboot, req->useStopEarly());
                vs->addPval(pval);
            }
        }

    //todo?
    //outputPvals(results, req.outputDir);

    return true;
}


class ParallelProcess
{
private:
    Request* req;
    SampleInfo* sampleInfo;

    std::vector<std::string> lines;
    std::vector<Variant> variants;
    std::vector<VariantSet*> pointers;

    VariantSet leftover;

    bool collapsing;
    bool parsing;
    bool testing;

    std::future<std::vector<Variant>> futureVariants;
    std::future<std::vector<VariantSet>> futureVariantSets;
    std::future<bool> testingDone;

public:

    ParallelProcess(Request *r, SampleInfo *si) : req(r), sampleInfo(si) {
        collapsing = false;
        parsing = false;
        testing = false;
    }

    inline bool isRunning(){ return collapsing || parsing || testing; }
    inline bool isParsing(){ return parsing; }
    inline bool isCollapsing(){ return collapsing; }
    inline bool isTesting(){ return testing; }

    //---------------------------------------------------
    void parseAndFilter(std::vector<std::string> & lns){
        parsing = true;
        this->lines = lns;

        futureVariants = std::async(std::launch::async,
                [this] { return constructVariants(req, sampleInfo, lines); }) ;
    }
    inline bool isParseDone(){
        return parsing && futureVariants.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }

    inline std::vector<Variant> getFilterResults(){
        parsing = false;
        lines.clear();
        return futureVariants.get();
    }
    //---------------------------------------------------

    void collapse(std::vector<Variant>& v, VariantSet& leftovers){
        collapsing = true;
        this->variants = v;
        this->leftover = leftovers;
        futureVariantSets = std::async(std::launch::async,
                [this] { return collapseVariants(req, variants, leftover); }) ;
    }
    inline bool isCollapseDone(){
        return collapsing && futureVariantSets.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }
    inline std::vector<VariantSet> getCollapseResults(){
        collapsing = false;
        variants.clear();
        return futureVariantSets.get();
    }
    //---------------------------------------------------

    void beginAssociationTest(std::vector<VariantSet*>& variants){
        testing = true;
        this->pointers = variants;

        testingDone = std::async(std::launch::async,
                [this] { return testBatch(req, sampleInfo, pointers); }) ;
    }
    inline bool isTestingDone(){
        return testing && testingDone.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }
    inline void setDone(){
        testing = false;
        pointers.clear();
    }
};

void processVCF(Request &req, SampleInfo &sampleInfo, std::vector<VariantSet>* results) {

    size_t nthreads = static_cast<size_t>(req.getNumberThreads());
    std::vector<ParallelProcess> threads;
    for(size_t i = 0; i < nthreads; i++)
        threads.emplace_back(&req, &sampleInfo);

    std::queue<ParallelProcess*> filterOrder;

    File vcf;
    vcf.open(req.getVCFDir());

    size_t totalLineCount = 0;
    size_t batchSize = static_cast<size_t>(req.getBatchSize());

    std::vector<std::string> lines;
    std::vector<Variant> constructedVariants;
    std::vector<VariantSet*> readyToRun;

    VariantSet leftover;
    bool readyToCollapse = false;

    //skips header
    extractHeaderLine(vcf);
    bool allParsingDone = false;
    bool allCollapsingDone = false;
    bool allTestingDone = false;

    while (!allParsingDone || !allCollapsingDone || !allTestingDone){

        //read in line if batch not full
        if(vcf.hasNext() && lines.size() < batchSize){
            lines.emplace_back(vcf.nextLine());
            totalLineCount++;

            if(totalLineCount % 100 == 0){
                if(totalLineCount == 0)
                    printInfo("Parsing VCF file...");
                else
                    printInfo(std::to_string(totalLineCount) + " variant lines have been parsed so far.");
            }
        }

        //check if parse thread is done
        if(filterOrder.size() > 0 && filterOrder.front()->isParseDone()){
            std::vector<Variant> v = filterOrder.front()->getFilterResults();
            filterOrder.pop();
            constructedVariants.insert(constructedVariants.end(), v.begin(), v.end());

            if(lines.size() == 0 && filterOrder.size() == 0 && !vcf.hasNext())
                allParsingDone = true;
        }

        bool collapseRunning = false;
        //check if collapse thread is done
        for(size_t m = 0; m < nthreads; m++){
            if(threads[m].isCollapseDone()){
                std::vector<VariantSet> v = threads[m].getCollapseResults();

                leftover = v.back();
                v.pop_back();

                results->insert(results->end(), v.begin(), v.end());

                size_t n = results->size() - 1;
                for(size_t i = 0; i < v.size(); i++)
                    readyToRun.emplace_back(&results->at(n-i));
                readyToCollapse = true;
            }
            collapseRunning = collapseRunning || threads[m].isCollapsing();
        }

        if(!collapseRunning && allParsingDone && constructedVariants.size() == 0){
            allCollapsingDone = true;
            if(leftover.size() > 0){
                results->push_back(leftover);
                readyToRun.push_back(&results->back());
            }
        }

        //check if testing thread is done
        bool allDone = true;
        for(size_t m = 0; m < nthreads; m++){

            if(threads[m].isTestingDone())
                threads[m].setDone();

            if(threads[m].isRunning())
                allDone = false;
        }

        if(allCollapsingDone && readyToRun.size() == 0 && allDone){
            allTestingDone = true;
            break;
        }

        //set up parse thread
        if(lines.size() >= batchSize || (!vcf.hasNext() && lines.size() > 0)){
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

        //set up collapse thread
        if(readyToCollapse && (constructedVariants.size() >= batchSize ||
                (constructedVariants.size() > 0 && allParsingDone))){
            for(size_t m = 0; m < nthreads; m++){
                 if(!threads[m].isRunning()){
                     threads[m].collapse(constructedVariants, leftover);
                     filterOrder.push(&threads[m]);
                     constructedVariants.clear();
                     readyToCollapse=false;

                     break;
                 }
             }
        }

        //set up test thread
        if(readyToRun.size() >= batchSize || (readyToRun.size() > 0 && allParsingDone)){
            for(size_t m = 0; m < nthreads; m++){
                 if(!threads[m].isRunning()){
                     threads[m].beginAssociationTest(readyToRun);
                     readyToRun.clear();
                     break;
                 }
             }
        }
    }

    printInfo("A total of " + std::to_string(totalLineCount) + " variants were parsed from the VCF file.");

    return;
}

