#include "../vikNGS.h"
#include "../Log.h"
#include "CLI11.h"

#include <iomanip>


std::string lower (std::string str){
    for(int i = 0; i < str.size(); i++)
        str[i] = std::tolower(str[i]);

    return str;
}


int main(int argc, char* argv[]) {

    CLI::App app{ "vikNGS Variant Association Toolkit" };

    std::string filename = "default";

    // -------------------------------------
    bool common = true;
    CLI::Option *c = app.add_flag("-c,--common", common, "Perform a common variant association test (default)");
    std::string rare = "";
    CLI::Option *r = app.add_option("-r,--rare", rare, "Perform a rare variant association test");

    std::string gt = "";
    CLI::Option *g = app.add_option("-g,--genotype", rare, "Genotype to use");

    // -------------------------------------

    // -------------------------------------
    int threads = 1;
    CLI::Option *t = app.add_option("-t,--threads", threads, "Number of threads", 1);
	t->check(CLI::Range(1, 2147483647));

    int nboot = 1;
    CLI::Option *n = app.add_option("-n,--boot", nboot, "Number of bootstrap iterations to calculate");
    n->check(CLI::Range(0, 2147483647));

    bool stopEarly = false;
    CLI::Option *s = app.add_flag("-s,--stop", stopEarly, "Stop bootstrapping if p-value looks to be > 0.05");
    // -------------------------------------

    // -------------------------------------
    std::string bedDir;
    CLI::Option *b = app.add_option("-b,--bed", bedDir, "Specify a directory of a BED file for collapsing variants");
    b->check(CLI::ExistingFile);

    bool byGene;
    CLI::Option *gene = app.add_option("--gene", byGene, "Collapse variants by gene if BED file specified (default)");
    bool byExon;
    CLI::Option *exon = app.add_option("--exon", byExon, "Collapse variants by exon if BED file specified");
    //bool byCoding;
    //CLI::Option *coding = app.add_option("--coding", byCoding, "Collapse variants by coding if BED file specified");
    int collapse = -1;
    CLI::Option *k = app.add_option("-k,--collapse", collapse, "Collapse every k variants (rare only)");
    k->check(CLI::Range(2, 2147483647));
    // -------------------------------------

    // -------------------------------------
    double maf = 0.05;
    CLI::Option *m = app.add_option("-m,--maf", maf, "Minor allele frequency cut-off (common-rare threshold)", 0.05);
    m->check(CLI::Range(0.0, 1.0));

    int depth = 30;
    CLI::Option *d = app.add_option("-d,--depth", depth, "Read depth cut-off (low-high read depth threshold)", 30);
    d->check(CLI::Range(1, 2147483647));

    double missing = 0.1;
    CLI::Option *x = app.add_option("-x,--missing", missing, "Missing data cut-off (variants with a proportion of missing data more than this threshold will not be tested)", 0.1);
    x->check(CLI::Range(0.0, 0.5));

    bool mustPass = false;
    CLI::Option *a = app.add_flag("-p,--pass", mustPass, "Include only variants which have PASS in the FILTER column");
	
    std::string chrFilter = "na";
    CLI::Option *chr = app.add_option("--chr", chrFilter, "Only include variants on this chromosome");

    int from = -1;
    CLI::Option *p1 = app.add_option("--from", from, "Only include variants with POS larger than this value");
    p1->check(CLI::Range(0, 2147483647));

    int to = -1;
    CLI::Option *p2 = app.add_option("--to", to, "Only include variants with POS smaller than this value");
    p2->check(CLI::Range(0, 2147483647));

    int batch = 1000;
    CLI::Option *h = app.add_option("-a,--batch", batch, "Processes VCF in batches of this many variants");
    h->check(CLI::Range(1, 2147483647));

    bool showFiltered = false;
    CLI::Option *filt = app.add_flag("--explain-filter", showFiltered, "Output explaination for filtered variants");

    // -------------------------------------

    // -------------------------------------
    std::string vcfDir;
    CLI::Option *v = app.add_option("vcf,-v,--vcf", vcfDir, "Specify a directory of a multisample VCF file (required)");
    v->required();
    v->check(CLI::ExistingFile);

    std::string sampleDir;
    CLI::Option *i = app.add_option("sample,-i,--sample", sampleDir, "Specify a directory of a TXT file containing sample information (required)");
    i->required();
    i->check(CLI::ExistingFile);

    std::string outputDir = ".";
    CLI::Option *o = app.add_option("-o,--out", outputDir, "Specify a directory for output (default = current directory)", ".");
    o->check(CLI::ExistingDirectory);

    // -------------------------------------

    CLI11_PARSE(app, argc, argv);

    printInfo("Starting vikNGS...");
    printInfo("VCF file: " + vcfDir);
    printInfo("Sample info file: " + sampleDir);
    printInfo("Output directory: " + outputDir);

    Request req = getDefaultRequest();
    req.setInputFiles(vcfDir, sampleDir);

    if(outputDir.back() == '/')
        outputDir.pop_back();

    req.setOutputDir(outputDir);

    printInfo("Minor allele frequency threshold: " + std::to_string(maf));
    req.setMafCutOff(maf);

    printInfo("Read depth threshold: " + std::to_string(depth));
    req.setHighLowCutOff(depth);

    printInfo("Missing data threshold: " + std::to_string(missing));
    req.setMissingThreshold(missing);

    if(chrFilter != "na"){
        printInfo("Analyzing variants on chromosome " + chrFilter);
        req.setChromosomeFilter(chrFilter);
    }
    if(from > 0){
        printInfo("Analyzing variants with POS greater than " + std::to_string(from));
        req.setMinPos(from);
    }
    if(from > 0){
        printInfo("Analyzing variants with POS less than " + std::to_string(to));
        req.setMaxPos(to);
    }

    if(mustPass)
        printInfo("Remove variants which do not PASS");

    req.setMustPASS(mustPass);

    if(nboot > 1){
        req.setBootstrap(nboot);
        req.setStopEarly(stopEarly);
        if(stopEarly)
            printInfo("Using " + std::to_string(nboot) + " bootstrap iterations and early stopping");
        else
            printInfo("Using " + std::to_string(nboot)  + " bootstrap iterations");
    }

    if(lower(rare) == "cast"){
        if(lower(gt) == "call"){
            req.addTest(TestSettings(GenotypeSource::CALL, Statistic::CAST, Variance::REGULAR));
            printInfo("Preparing to run rare variant association (CAST p-values, called genotypes)...");
        }
        else if(lower(gt) == "vcf"){
            req.addTest(TestSettings(GenotypeSource::VCF_CALL, Statistic::CAST, Variance::RVS));
            printInfo("Preparing to run rare variant association (CAST p-values, VCF GT)...");
        }
        else{
            req.addTest(TestSettings(GenotypeSource::EXPECTED, Statistic::CAST, Variance::REGULAR));
            printInfo("Preparing to run rare variant association (CAST p-values, expected GT/vRVS)...");
        }
    }
    else if(lower(rare) == "skat"){
        if(lower(gt) == "call"){
            req.addTest(TestSettings(GenotypeSource::CALL, Statistic::SKAT, Variance::REGULAR));
            printInfo("Preparing to run rare variant association (SKAT p-values, called genotypes)...");
        }
        else if(lower(gt) == "vcf"){
            req.addTest(TestSettings(GenotypeSource::VCF_CALL, Statistic::SKAT, Variance::RVS));
            printInfo("Preparing to run rare variant association (SKAT p-values, VCF GT)...");
        }
        else{
            req.addTest(TestSettings(GenotypeSource::EXPECTED, Statistic::SKAT, Variance::REGULAR));
            printInfo("Preparing to run rare variant association (SKAT p-values, expected GT/vRVS)...");
        }
    }
    else{
        if(lower(gt) == "call"){
            req.addTest(TestSettings(GenotypeSource::CALL, Statistic::COMMON, Variance::REGULAR));
            printInfo("Preparing to run common variant association (called genotypes)...");
        }
        else if(lower(gt) == "vcf"){
            req.addTest(TestSettings(GenotypeSource::VCF_CALL, Statistic::SKAT, Variance::RVS));
            printInfo("Preparing to run common variant association (VCF GT)...");
        }
        else{
            req.addTest(TestSettings(GenotypeSource::EXPECTED, Statistic::SKAT, Variance::REGULAR));
            printInfo("Preparing to run common variant association (expected GT/vRVS)...");
        }
    }

    if(collapse > 2){
        printInfo("Collapse every " + std::to_string(collapse) + " variants");
        req.setCollapse(collapse);
    }

    if(bedDir.size() > 0){

        req.setCollapseFile(bedDir);
        printInfo("BED file: " + bedDir);

        if(byGene){
            printInfo("Collapse variants along genes");
            req.setCollapseGene();
        }
        else if(byExon){
            printInfo("Collapse variants along exons");
            req.setCollapseExon();
        }
    }

    printInfo("Batch size: " + std::to_string(batch));
    req.setBatchSize(batch);

    if(threads != 1){
        printInfo("Using " + std::to_string(threads) + " threads");
    }
    req.setNumberThreads(threads);

    req.setKeepFiltered(showFiltered);

    startVikNGS(req);
    return 0;
}

