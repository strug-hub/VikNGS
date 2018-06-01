#include "RVS.h"

#include <iostream>  
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>


void generateForR(MatrixXd X, VectorXd Y, MatrixXd Z, VectorXd G, MatrixXd P, std::map<int, int> readGroup) {
	std::ofstream Xfile("C:/Users/Scott/Desktop/RVS-master/example/X.txt");
	std::ofstream Yfile("C:/Users/Scott/Desktop/RVS-master/example/Y.txt");
	std::ofstream Pfile("C:/Users/Scott/Desktop/RVS-master/example/P.txt");
	std::ofstream Mfile("C:/Users/Scott/Desktop/RVS-master/example/M.txt");
	std::ofstream Zfile("C:/Users/Scott/Desktop/RVS-master/example/Z.txt");

	int precise = 18;

	if (Mfile.is_open())
	{
		Mfile << "ID\thrg\n";

		for (size_t i = 0; i < G.rows(); i++) {

			Mfile << G[i];
			Mfile << '\t';
			Mfile << readGroup[G[i]];
			Mfile << '\n';
		}
		Mfile.close();
	}


	if (Yfile.is_open())
	{
		for (size_t i = 0; i < Y.rows(); i++) {
			Yfile << std::setprecision(precise) << Y[i];
			Yfile << '\n';
		}
		Yfile.close();
	}

	if (Xfile.is_open())
	{
		for (size_t i = 0; i < X.cols(); i++) {
			for (size_t j = 0; j < X.rows(); j++) {

				if (std::isnan(X(j, i)))
					Xfile << "NA";
				else
					Xfile << std::setprecision(precise) << X(j, i);

				if (j < X.rows() - 1)
					Xfile << '\t';
			}
			Xfile << '\n';
		}

		Xfile.close();
	}

	if (Pfile.is_open())
	{
		for (size_t i = 0; i < P.rows(); i++) {
			Pfile << std::setprecision(precise) << P(i, 0);
			Pfile << '\t';
			Pfile << std::setprecision(precise) << P(i, 1);
			Pfile << '\t';
			Pfile << std::setprecision(precise) << P(i, 2);
			Pfile << '\n';
		}

		Pfile.close();
	}

	if (Zfile.is_open())
	{
		for (size_t i = 0; i < Z.rows(); i++) {
			for (size_t j = 0; j < Z.cols(); j++) {
				Zfile << std::setprecision(precise) << Z(i, j);

				if (j < Z.cols() - 1)
					Zfile << '\t';
			}

			Zfile << '\n';
		}

		Zfile.close();
	}
}


void run(Request req) {


  createFile(req.outputDir);

  printInfo("Parsing files...");
  TestInput input = parseAndFilter(req);

  if(input.hasCovariates())
    printInfo(std::to_string(input.countCovariates()) + " covariates parsed");

  if (req.useCommon()) {
    
    printInfo("Starting tests...");
    std::vector<Variant> variants = runCommonTest(req, input);
    
    printInfo("Common Test p-values");
    for (size_t i = 0; i < variants.size(); i++) {
      std::cout << variants[i].toString() << "\t" << variants[i].pvalue;
      std::cout << '\n';
    }
    
    outputPvals(variants, req.outputDir);
  }
  
  if (!req.useCommon()) {
    
    printInfo("Starting tests...");
    std::vector<Variant> variants = runRareTest(req, input);
    
    printInfo("Rare Test p-values");
    for (size_t i = 0; i < variants.size(); i++) {
      std::cout << variants[i].toString() << "\t" << variants[i].pvalue;
      std::cout << '\n';
    }

    outputPvals(variants, req.outputDir, 5);
   }
}

#include "CLI11.hpp"

int main(int argc, char* argv[]) {

	CLI::App app{ "vikNGS Variant Association Toolkit" };

	std::string filename = "default";
	
	// -------------------------------------
	bool common = true;
	CLI::Option *c = app.add_flag("-c,--common", common, "Perform a common variant association test (default)");
	std::string rare = "";
	CLI::Option *r = app.add_option("-r,--rare", rare, "Perform a rare variant association test");

	// -------------------------------------

	// -------------------------------------
	int nthreads = 1;
	CLI::Option *t = app.add_option("-t,--threads", nthreads, "Number of threads", 1);
	t->check(CLI::Range(1, 2147483647));

	int nboot = 0;
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
	int collapse;
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

	bool all = false;
	CLI::Option *a = app.add_flag("-a,--all", all, "Include variants which do not have PASS in the FILTER column");
	
	int from = -1;
	CLI::Option *p1 = app.add_option("--from", from, "Only include variants with POS larger than this value");
	p1->check(CLI::Range(0, 2147483647));
	
	int to = -1;
	CLI::Option *p2 = app.add_option("--to", to, "Only include variants with POS smaller than this value");
	p2->check(CLI::Range(0, 2147483647));

	int batch = 1000;
	CLI::Option *h = app.add_option("-h,--batch", batch, "Processes VCF in batches of this many variants");
	h->check(CLI::Range(1, 2147483647));
	
	// -------------------------------------

	// -------------------------------------
	std::string vcfDir;
	CLI::Option *i = app.add_option("vcf,-i,--vcf", vcfDir, "Specify a directory of a multisample VCF file (required)");
	i->required();
	i->check(CLI::ExistingFile);

	std::string sampleDir;
	CLI::Option *g = app.add_option("sample,-g,--sample", sampleDir, "Specify a directory of a TXT file containing sample information (required)");
	g->required();
	g->check(CLI::ExistingFile);

	std::string outputDir = ".";
	CLI::Option *o = app.add_option("-o,--out", outputDir, "Specify a directory for output (default = current directory)", ".");
	o->check(CLI::ExistingDirectory);

	// -------------------------------------

	CLI11_PARSE(app, argc, argv);

	printInfo("Starting vikNGS...");
	printInfo("Sample info file: " + sampleDir);
	printInfo("Output directory: " + outputDir);

	initializeRequest(vcfDir, sampleDir);
	setOutputDir(outputDir);

	setMAFCutoff(maf);
	printInfo("Minor allele frequency threshold: " + std::to_string(maf));

	setHighLowCutOff(depth);
	printInfo("Read depth threshold: " + std::to_string(depth));

	setMissingThreshold(missing);
	printInfo("Missing data threshold: " + std::to_string(missing));

	setMustPASS(all);
	if(a->count() > 0)
		printInfo("Retain variants which do not PASS");
	else
		printInfo("Remove variants which do not PASS");
	
	  if(from != -1 && to != -1){
	    setMinPos(from);
	    setMaxPos(to);
	    printInfo("Analyzing variants between POS " + std::to_string(from) + " and " + std::to_string(to));
	  }
	  else if(from != -1){
	    setMinPos(from);
	    printInfo("Analyzing variants with POS greater than " + std::to_string(from));
	  }
	  else if(to != -1){
	    setMinPos(to);
	    printInfo("Analyzing variants with POS less than " + std::to_string(to));
	  }

	if(r->count() > 0){

		if(rare=="calpha"){
		  useRareTest(rare);
		  printInfo("Preparing to run rare variant association (C-alpha p-values)...");
		}

		else if(rare=="cast"){
		  useRareTest(rare);
		  printInfo("Preparing to run rare variant association (CAST p-values)...");
		}
		else{
		  useRareTest("cast");
		  printInfo("Preparing to run rare variant association (defaulting to CAST p-values)...");
		}

		if(k->count() > 0 && b->count() < 1){
		  setCollapse(collapse);
		  printInfo("Collapsing every " + std::to_string(collapse) + " variants");
		}

	}
	else{
		useCommonTest();
		printInfo("Preparing to run common variant association...");
	}

	if(nboot > 1){
		useBootstrap(nboot);

		setStopEarly(stopEarly);
		if(stopEarly)
			printInfo("Using " + std::to_string(nboot) + " bootstrap iterations and early stopping");
		else
			printInfo("Using " + std::to_string(nboot) + " bootstrap iterations");

	}

	setNumberThreads(nthreads);
	if(nthreads != 1){
		printInfo("Using " + std::to_string(nthreads) + " threads");
	}

	if(b->count() > 0){

		setCollapseFile(bedDir);
		printInfo("BED file: " + bedDir);

		if(exon->count() > 0){
			setCollapseExon();
			printInfo("Collapse variants along exons");
		}/*
		else if(coding->count() > 0){
			setCollapseCoding();
			printInfo("Collapse variants along coding regions");
		}*/
		else{
			setCollapseGene();
			printInfo("Collapse variants along genes");
		}
	}

	if(h->count() > 0)
       		setBatchSize(batch);


	Request req = getRequest();
	run(req);

	return 0;
}

