#include "Parser.h"
#include "Filter.h"
#include "../Test/Test.h"
#include "File.h"
#include "../Output/OutputHandler.h"
#include "../Request.h"
#include "../SampleInfo.h"
#include "../vikNGS.h"
#include "../Log.h"

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

        if(STOP_RUNNING_THREAD)
            return variants;

        //extract to the FILTER column
        std::vector<std::string> info = splitString(lines[i], VCF_SEP, FILTER + 1);

        if(info.size() < FORMAT)
            continue;

        Filter filter = filterByVariantInfo(req, info[CHROM], info[POS], info[REF], info[ALT], info[FILTER]);

        if(filter == Filter::IGNORE)
            continue;

        Variant variant;

        if(filter == Filter::VALID){
            std::vector<std::string> columns = splitString(lines[i], VCF_SEP);
            variant = constructVariant(columns, calculateExpected, calculateCalls, getVCFCalls);

            if(variant.isValid()){
                VectorXd Y = sampleInfo->getY();
                filter = filterByGenotypes(req, variant, Y, sampleInfo->getFamily());
            }
            else
                continue;

        }
        else
            try{ variant = Variant(info[CHROM], std::stoi(info[POS]), info[ID], info[REF], info[ALT]); }
            catch(...){ continue; }

        variant.setFilter(filter);
        variants.push_back(variant);
    }
    variants.shrink_to_fit();
    return variants;
}

int findInterval(IntervalSet * is, std::string chr, int pos, int searchHint){
    size_t hint = static_cast<size_t>(searchHint);
    std::vector<Interval>* intervals = is->get(chr);
    if(intervals->size() < 1) {
        printBedIDWarning(chr);
        return -1;
    }
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
        if(left == right)
            break;
    }

    if (intervals->at(middle).isIn(pos))
          return static_cast<int>(middle);

    return -1;
}

std::vector<VariantSet> collapseVariants(Request* req, std::vector<Variant>& variants, VariantSet& leftover){

    CollapseType collapse = req->getCollapseType();
    int k = req->getCollapseSize();
    std::vector<VariantSet> variantSet;

    if(collapse == CollapseType::NONE){
        if(leftover.size() > 0)
            variantSet.push_back(leftover);
        for(size_t i = 0; i < variants.size(); i++)
            variantSet.emplace_back(VariantSet(variants[i]));
    }

    else if(collapse == CollapseType::COLLAPSE_K){

        size_t startFrom = 0;
        while(leftover.validSize() < k && startFrom < variants.size()){
            leftover.addVariant(variants[startFrom]);
            startFrom++;
        }
        variantSet.push_back(leftover);

        for(size_t i = startFrom; i < variants.size(); i++){
            if(variantSet.back().validSize() == k)
                variantSet.emplace_back(VariantSet(variants[i]));
            else
                variantSet.back().addVariant(variants[i]);
        }
    }
    else if(collapse == CollapseType::COLLAPSE_EXON || collapse == CollapseType::COLLAPSE_GENE){

        IntervalSet* is = req->getIntervals();

        size_t startFrom = 0;
        int hint = -1;

        while(startFrom < variants.size() && leftover.size() < 1){
            int index = findInterval(is, variants[startFrom].getChromosome(), variants[startFrom].getPosition(), -1);
            if(index >= 0){
                leftover.addVariant(variants[startFrom]);
                hint = index + 1;
                leftover.setInterval(&is->get(variants[startFrom].getChromosome())->at(index));
            }
            startFrom++;
        }

        if(leftover.size() > 0)
            variantSet.push_back(leftover);

        for(size_t i = startFrom; i < variants.size(); i++){

            if(variantSet.back().isIn(variants[i]) && variantSet.back().size() < k)
                variantSet.back().addVariant(variants[i]);
            else{
                int index = findInterval(is, variants[i].getChromosome(), variants[i].getPosition(), hint);
                if(index < 0)
                    continue;

                VariantSet newSet(variants[i]);
                newSet.setInterval(&is->get(variants[i].getChromosome())->at(index));
                variantSet.push_back(newSet);
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

    for(int i = 0; i < variants.size(); i++){
        for(Test t : req->getTests()){
            if(variants[i]->validSize() > 0){
                double pval = runTest(sampleInfo, variants[i], t, nboot, req->useStopEarly());
                variants[i]->addPval(pval);
            }
        }
    }

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
    void parseAndFilter(std::deque<std::string>& lns, size_t n){
        n = (lns.size() < n) ? lns.size() : n;
        if(n < 1) return;
        lines.clear(); lines.reserve(n);

        parsing = true;

        for(size_t i = 0; i < n; i++){
            lines.push_back(lns.front());
            lns.pop_front();
        }

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

    void collapse(std::deque<Variant>& v, VariantSet& leftovers, size_t n){
        n = (v.size() < n) ? v.size() : n;
        if(n < 1) return;
        variants.clear(); variants.reserve(n);

        collapsing = true;

        for(size_t i = 0; i < n; i++){
            variants.push_back(v.front());
            v.pop_front();
        }

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

    void beginAssociationTest(std::deque<VariantSet*>& vs, size_t n){
        n = (vs.size() < n) ? vs.size() : n;
        if(n < 1) return;
        pointers.clear(); pointers.reserve(n);

        testing = true;

        for(size_t i = 0; i < n; i++){
            pointers.push_back(vs.front());
            vs.pop_front();
        }

        testingDone = std::async(std::launch::async,
                [this] { return testBatch(req, sampleInfo, pointers); }) ;
    }
    inline bool isTestingDone(){
        return testing && testingDone.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }

    inline void setDone(){
        parsing = false;
        collapsing = false;
        testing = false;
        pointers.clear();
        lines.clear();
        variants.clear();
        return;
    }

};

std::vector<VariantSet> processVCF(Request &req, SampleInfo &sampleInfo, size_t& totalLineCount) {
    std::vector<VariantSet> results;

    size_t nthreads = std::max(1, req.getNumberThreads());
    std::vector<ParallelProcess> threads;
    for(size_t i = 0; i < nthreads; i++)
        threads.emplace_back(&req, &sampleInfo);

    std::queue<ParallelProcess*> parseOrder;
    std::queue<ParallelProcess*> collapseOrder;

    File vcf;
    vcf.open(req.getVCFDir());

    totalLineCount = 0;
    size_t batchSize = static_cast<size_t>(req.getBatchSize());

    std::deque<std::string> lines;
    std::deque<Variant> constructedVariants;
    std::deque<VariantSet> collapsedVariants;
    std::deque<VariantSet*> readyToRun;

    VariantSet leftover;

    //skips header
    extractHeaderLine(vcf);
    bool allParsingDone = false;
    bool allCollapsingDone = false;
    bool allTestingDone = false;

    while (!allParsingDone || !allCollapsingDone || !allTestingDone){

        if(STOP_RUNNING_THREAD)
            break;

        //read in line if batch not full
        if(vcf.hasNext() && lines.size() < batchSize){
            lines.emplace_back(vcf.nextLine());
            totalLineCount++;

            if(totalLineCount % batchSize == 0){
                if(totalLineCount == 0)
                    printInfo("Parsing VCF file...");
                else
                    printInfo(std::to_string(totalLineCount) + " variant lines have been parsed so far.");
            }
        }

        //check if parse thread is done
        if(parseOrder.size() > 0 && parseOrder.front()->isParseDone()){
            std::vector<Variant> v = parseOrder.front()->getFilterResults();
            parseOrder.pop();

            std::vector<Variant> filtered;
            for(size_t i = 0; i < v.size(); i++){
                if(v[i].isValid())
                    constructedVariants.push_back(v[i]);
                else if(req.shouldKeepFiltered())
                    filtered.push_back(v[i]);
            }

            if(filtered.size() > 0)
                outputFiltered(filtered, req.getOutputDir());

        }
        if(lines.size() == 0 && parseOrder.size() == 0 && !vcf.hasNext()){
            if(!allParsingDone)
                printInfo("A total of " + std::to_string(totalLineCount) + " variants were parsed from the VCF file.");
            allParsingDone = true;
        }

        //check if collapse thread is done
        if(collapseOrder.size() > 0 && collapseOrder.front()->isCollapseDone()){
            std::vector<VariantSet> v = collapseOrder.front()->getCollapseResults();
            collapseOrder.pop();

            if(v.size() > 0){
                leftover = v.back();
                v.pop_back();

                for(size_t i = 0; i < v.size(); i++){
                    collapsedVariants.push_back(v[i]);
                    readyToRun.push_back(&collapsedVariants.back());
                }
            }
       }

        if(collapseOrder.size() == 0 && allParsingDone && constructedVariants.size() == 0){
            allCollapsingDone = true;
            if(leftover.size() > 0){
                collapsedVariants.push_back(leftover);
                readyToRun.push_back(&collapsedVariants.back());
                VariantSet empty; leftover = empty;
            }
        }

        //check if testing thread is done
        bool allDone = true;
        for(size_t m = 0; m < nthreads; m++){

            if(threads[m].isTestingDone()){
                threads[m].setDone();
                while(collapsedVariants.size() > 0 && collapsedVariants.front().nPvals() > 0 &&
                      (collapsedVariants.size() < batchSize || collapsedVariants[batchSize-1].nPvals() > 0)){
                    results.push_back(collapsedVariants.front());
                    collapsedVariants.pop_front();
                    if(!req.shouldRetainGenotypes())
                        results.back().shrink();
                }
            }

            if(threads[m].isRunning())
                allDone = false;
        }

        if(allCollapsingDone && collapsedVariants.size() == 0 && allDone){
            allTestingDone = true;
            break;
        }

        //set up test thread
        if(readyToRun.size() >= batchSize || (readyToRun.size() > 0 && allCollapsingDone)){
            for(size_t m = 0; m < nthreads; m++){
                 if(!threads[m].isRunning()){
                     int size = batchSize;
                     if(allCollapsingDone) size = readyToRun.size();
                     threads[m].beginAssociationTest(readyToRun, size);
                     break;
                 }
             }
        }

        //set up collapse thread
        if(collapseOrder.size() == 0 && (constructedVariants.size() >= batchSize ||
                (constructedVariants.size() > 0 && allParsingDone))){
            for(size_t m = 0; m < nthreads; m++){
                 if(!threads[m].isRunning()){
                     threads[m].collapse(constructedVariants, leftover, batchSize);
                     collapseOrder.push(&threads[m]);
                     VariantSet empty; leftover = empty;
                     break;
                 }
             }
        }


        //set up parse thread
        if(lines.size() >= batchSize || (!vcf.hasNext() && lines.size() > 0)){
            //start new thread
            for(size_t m = 0; m < nthreads; m++){
                 if(!threads[m].isRunning()){
                     threads[m].parseAndFilter(lines, batchSize);
                     parseOrder.push(&threads[m]);
                     break;
                 }
             }
        }

        //set up small test thread
        if(readyToRun.size() > 0){
            if(req.getCollapseType() == CollapseType::COLLAPSE_EXON ||
                    req.getCollapseType() == CollapseType::COLLAPSE_GENE){
                for(size_t m = 0; m < nthreads; m++){
                     if(!threads[m].isRunning()){
                         int size = (readyToRun.size() < 3) ? readyToRun.size() : 3;
                         threads[m].beginAssociationTest(readyToRun, size);
                         break;
                     }
                }
            }
        }

    }

    while(true){
        bool threadsDone = true;
        for(size_t m = 0; m < nthreads; m++){
             if(threads[m].isRunning()){
                 threadsDone = false;
                 if(threads[m].isTestingDone())
                     threads[m].setDone();
                 if(threads[m].isCollapseDone())
                     threads[m].getCollapseResults();
                 if(threads[m].isParseDone())
                     threads[m].getFilterResults();
             }
        }

        if(threadsDone)
            break;
    }

    return results;
}
