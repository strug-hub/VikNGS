#include "stdafx.h"
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

				if (isnan(X(j, i)))
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

SimulationRequestGroup newSimulationRequestGroup(int groupID,
	std::string n, std::string caseControl, std::string highLow, 
	std::string meanDepth, std::string sdDepth) {

	std::string id = std::to_string(groupID);

	SimulationRequestGroup requestGroup;

	try { requestGroup.n = std::abs(std::stoi(n)); }
	catch (...) { throw std::invalid_argument("group " + id + " sample size,integer"); }
	if (caseControl != "case" && caseControl != "control")
		throw std::invalid_argument("group " + id + " cohort,one of 'case' or 'control'");
	requestGroup.isCase = (caseControl == "case"); 
	if (highLow != "high" && highLow != "low")
		throw std::invalid_argument("group " + id + " read depth,one of 'high' or 'low'");
	requestGroup.isHrg = (highLow == "high"); 

	try { requestGroup.meanDepth = std::abs(std::stod(meanDepth)); }
	catch (...) { throw std::invalid_argument("group " + id + " mean read depth,decimal"); }
	try { requestGroup.sdDepth = std::abs(std::stod(sdDepth)); }
	catch (...) { throw std::invalid_argument("group " + id + " read depth standard deviation,decimal"); }
	
	return requestGroup;
}

void validateSimulationRequest(SimulationRequest request) {

	int nsample = 0;
	for (SimulationRequestGroup group : request.groups)
		nsample += group.n;

	if (nsample > request.npop)
		throw std::domain_error("Population size (value given: " +
			std::to_string(request.npop) +
			") should be greater than sample size (summed across groups: " +
			std::to_string(nsample) + ").");

	if (request.npop <= 0)
		throw std::domain_error("Population size (value given: " +
			std::to_string(request.npop) +
			") should be greater than zero.");

	if (request.nsnp <= 0)
		throw std::domain_error("Number of variants (value given: " +
			std::to_string(request.nsnp) +
			") should be greater than zero.");

	if (request.prevalence > 1)
		throw std::domain_error("Revelance (value given: " +
			std::to_string(request.prevalence) +
			") should be a value between 0 and 1.");

	if (request.maf > 0.5)
		throw std::domain_error("Minor allele frequency (value given: " +
			std::to_string(request.maf) +
			") should be a value between 0 and 0.5.");

	if (request.groups.size() < 2)
		throw std::domain_error("Simulation requires at least two groups (" +
			std::to_string(request.groups.size()) + " given).");


	//if (!request.useCommonTest && request.nsnp < 5)
		//throw std::domain_error("Rare test requires number of variants to be at least 5 (value given: " +
			//std::to_string(request.nsnp) + ").");

	if (request.useBootstrap && request.nboot < 1)
		throw std::domain_error("Nuber of bootstrap iterations should be at least 1 (value given: " +
			std::to_string(request.nboot) + ").");

	int ncase = 0;
	int ncontrol = 0;
	for (SimulationRequestGroup g : request.groups) {
		if (g.isCase)
			ncase++;
		if (!g.isCase)
			ncontrol++;

		if(g.n < 1)
			throw std::domain_error("Each group should have a sample size of at least 1.");
	}

	if (ncase < 1 || ncontrol < 1) 		
		throw std::domain_error("Case/control simulation requires at least one case and one control group (" +
			std::to_string(ncase) + " case and " +
			std::to_string(ncontrol) + " control groups given).");

}


SimulationRequest newSimulationRequest(std::string npop, std::string prevalence,
	std::string nsnp, std::string me, std::string sde, std::string oddsRatio,
	std::string maf, std::vector<SimulationRequestGroup> groups,
	std::string test, bool useBootstrap, std::string nboot) {

	SimulationRequest request;

	try { request.npop = std::abs(std::stoi(npop)); }
	catch (...) { throw std::invalid_argument("population size,integer"); }
	try { request.prevalence = std::abs(std::stod(prevalence)); }
	catch (...) { throw std::invalid_argument("prevalence,decimal"); }
	try { request.nsnp = std::abs(std::stoi(nsnp)); }
	catch (...) { throw std::invalid_argument("number of variants,integer"); }
	try { request.me = std::abs(std::stod(me)); }
	catch (...) { throw std::invalid_argument("sequencing error rate,decimal"); }
	try { request.sde = std::abs(std::stod(sde)); }
	catch (...) { throw std::invalid_argument("sequencing error standard deviation,integer"); }	
	try { request.oddsRatio = std::abs(std::stod(oddsRatio)); }
	catch (...) { throw std::invalid_argument("odds ratio,decimal"); }
	try { request.maf = std::abs(std::stod(maf)); }
	catch (...) { throw std::invalid_argument("Minor allele frequency,decimal"); }

	request.test = test;
	request.useBootstrap = useBootstrap;

	if (useBootstrap) {
		try { request.nboot = std::abs(std::stoi(nboot)); }
		catch (...) { throw std::invalid_argument("# bootstrap iterations,integer"); }
	}

	request.groups = groups;

	try { validateSimulationRequest(request); }
	catch (...) { throw; }

	return request;
}

std::vector<double> startVikNGS(Request req) {
	
	std::vector<double> p;

	VectorXd Y, G; MatrixXd X, Z, P;
	std::map<int, int> readGroup;
	std::vector<std::vector<int>> interval;

	bool valid = parseAndFilter(req, X, Y, Z, G, readGroup, P, interval);

	bool useCovariates = Z.rows() > 0;

	std::vector<double> pval;


	if (!valid) {
		//TODO 
		return pval;
	}

    if (req.useCommon()) {

		if (req.useBootstrap && useCovariates)
			pval = runCommonTest(X, Y, Z, G, readGroup, P, req.nboot);
		else if(req.useBootstrap && !useCovariates)
			pval = runCommonTest(X, Y, G, readGroup, P, req.nboot);
		else if (!req.useBootstrap && useCovariates)
			pval = runCommonTest(X, Y, Z, G, readGroup, P);
		else
			pval = runCommonTest(X, Y, G, readGroup, P);

	}
	else {
		if (useCovariates)
			pval = runRareTest(X, Y, Z, G, readGroup, P, req.nboot);
		else
			pval = runRareTest(X, Y, G, readGroup, P, req.nboot);
	}

	std::string outputDir = "C:/Users/Scott/Desktop/out.txt";
	outputPvals(pval, outputDir);

	return pval;
}

std::vector<double> startSimulation(SimulationRequest req) {

	VectorXd Y, G; MatrixXd X, P;
	std::map<int, int> readGroup;

	simulate(req, X, Y, G, readGroup, P);

	std::vector<double> pval;

	if (req.test == "common"){
		
		if(req.useBootstrap)
			pval = runCommonTest(X, Y, G, readGroup, P, req.nboot);
		else
			pval = runCommonTest(X, Y, G, readGroup, P);

	}
	else
		pval = runRareTest(X, Y, G, readGroup, P, req.nboot, req.test);

    std::string outputDir = "C:/Users/Scott/Desktop/out.txt";
    outputPvals(pval, outputDir);

    return pval;
}

SimulationRequest testSimulationRequest() {
	std::vector<SimulationRequestGroup> groups;
	groups.push_back(newSimulationRequestGroup(0, "300", "control", "low", "10", "5"));
	groups.push_back(newSimulationRequestGroup(0, "200", "case", "high", "40", "7"));
	groups.push_back(newSimulationRequestGroup(0, "100", "control", "high", "50", "6"));

    return newSimulationRequest("3000", "0.1", "100", "0.01", "0.025", "1.4", "0.1", groups, "common", false, 0);
}



int main() {

    //TODO: take as input from command line
    //---------------------------------------
	bool simulation = false;

	///input files
    //std::string vcfDir = "C:/Users/Scott/Desktop/vcf/chr7_case_control.vcf";
    //std::string infoDir = "C:/Users/Scott/Desktop/vcf/sampleInfo.txt";

	std::string vcfDir = "C:/Users/Scott/Desktop/vcf/example_1000snps.vcf";
	std::string infoDir = "C:/Users/Scott/Desktop/vcf/sampleInfo.txt";

    std::string bedDir = "";
    //std::string bedDir = "C:/Users/Scott/Desktop/RVS-master/example/chr11.bed";

    std::string outputDir = "C:/Users/Scott/Desktop/out.txt";

	initializeRequest(vcfDir, infoDir);
	useCommonTest();
	setHighLowCutOff(30);
	
	//setCollapseFile(bedDir);
	//setCollapseCoding();

	setMAFCutoff(0.05);
	setMissingThreshold(0.1);
	setMustPASS(true);
	setOnlySNPs(true);
	Request req = getRequest();


	bool rvs = true;

    //---------------------------------------

    //TODO: check to see if file can be opened when another application is using it (excel)
    //TODO: test windows vs unix EOF characters, doesn't seem to work well with windows

    VectorXd Y, G; MatrixXd X, Z, P;
    std::map<int, int> readGroup;
    std::vector<std::vector<int>> interval;

    if (simulation) {
        simulate(testSimulationRequest(), X, Y, G, readGroup, P);
        std::vector<double> pvals = runCommonTest(X, Y, G, readGroup, P);

        std::cout << "Common Test p-values\n";
        for (size_t i = 0; i < pvals.size(); i++) {
            std::cout << pvals[i];
            std::cout << '\n';
        }

        outputPvals(pvals, outputDir);
    }
    else {
        bool valid = parseAndFilter(req, X, Y, Z, G, readGroup, P, interval);
        if (!valid) {
            return 0;
        }

        generateForR(X, Y, Z, G, P, readGroup);

        if (req.useCommon()) {
            std::vector<double> pvals = runCommonTest(X, Y, G, readGroup, P);
            //std::vector<double> pvals = runCommonTest(X, Y, Z, G, readGroup, P);
            //std::vector<double> pvals = runCommonTest(X, Y, Z, G, readGroup, P, 1000, true);

            std::cout << "Common Test p-values\n";
            for (size_t i = 0; i < pvals.size(); i++) {
                std::cout << pvals[i];
                std::cout << '\n';
            }

			outputPvals(pvals, outputDir);

        }
        else {

            std::vector<double> pval = runRareTest(X, Y, G, readGroup, P, 10);
            //std::vector<double> pval = runRareTest(X, Y, G, readGroup, P, Z, 20000, true);

            std::cout << "Rare Test p-values\n";
			for (size_t i = 0; i < pval.size(); i++) {
				std::cout << pval[i];
				std::cout << '\n';
			}

			outputPvals(pval, outputDir);

        }
    }



    //keep console open while debugging
    //TODO: be sure to remove eventually!
    std::cout << "\ndone...>";
    while (true) {}
    return 0;
}



