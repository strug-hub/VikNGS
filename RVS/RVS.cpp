#include "stdafx.h"
#include "RVS.h"

#include <iostream>  
#include <string>
#include <vector>

#include <iostream>
#include <fstream>
#include <iomanip>


void generateForR(MatrixXd X, VectorXd Y, MatrixXd Z, VectorXd G, MatrixXd P, std::map<int, int> readGroup) {
	std::ofstream Xfile("C:/Users/Scott/Desktop/RVS-master/example/X.txt");
	std::ofstream Yfile("C:/Users/Scott/Desktop/RVS-master/example/Y.txt");
	std::ofstream Pfile("C:/Users/Scott/Desktop/RVS-master/example/P.txt");
	std::ofstream Mfile("C:/Users/Scott/Desktop/RVS-master/example/M.txt");
	std::ofstream Zfile("C:/Users/Scott/Desktop/RVS-master/example/Z.txt");

	int precise = 12;

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

	if (request.lowerMAF > 0.5)
		throw std::domain_error("MAF lower bound (value given: " +
			std::to_string(request.lowerMAF) +
			") should be a value between 0 and 0.5.");
	if (request.upperMAF > 0.5)
		throw std::domain_error("MAF upper bound (value given: " +
			std::to_string(request.upperMAF) +
			") should be a value between 0 and 0.5.");

	if (request.lowerMAF > request.upperMAF)
		throw std::domain_error("MAF lower bound (value given: " +
			std::to_string(request.lowerMAF) +
			") should be less than upper bound (value given: " +
			std::to_string(request.upperMAF) + ")");

	if (request.groups.size() < 2)
		throw std::domain_error("Simulation requires at least two groups (" +
			std::to_string(request.groups.size()) + " given).");

	if (!request.useCommonTest && request.nsnp < 5)
		throw std::domain_error("Rare test requires number of variants to be at least 5 (value given: " +
			std::to_string(request.nsnp) + ").");

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

void validateRequest(Request request) {
	//TODO
}

Request newRequest(std::string vcfDir, std::string sampleDir, std::string bedDir,
	std::string highLowCutOff, bool collapseCoding, bool collapseExon,
	std::string mafCutoff, std::string missingThreshold, bool onlySNPs, bool mustPASS,
	bool useCommonTest, bool useBootstrap, std::string nboot) {

	Request request;

	request.vcfDir = vcfDir;
	request.sampleDir = sampleDir;
	request.bedDir = bedDir;
	request.collapseCoding = collapseCoding;
	request.collapseExon = collapseExon;
	request.onlySNPs = onlySNPs;
	request.mustPASS = mustPASS;
	request.useCommonTest = useCommonTest;
	request.useBootstrap = useBootstrap;


	try { request.highLowCutOff = std::abs(std::stoi(highLowCutOff)); }
	catch (...) { throw std::invalid_argument("read depth threshold,integer"); }
	try { request.mafCutoff = std::abs(std::stod(mafCutoff)); }
	catch (...) { throw std::invalid_argument("MAF cutoff,decimal"); }
	try { request.missingThreshold = std::abs(std::stod(missingThreshold)); }
	catch (...) { throw std::invalid_argument("missing threshold,decimal"); }

	if (useBootstrap) {
		try { request.nboot = std::abs(std::stoi(nboot)); }
		catch (...) { throw std::invalid_argument("# bootstrap iterations,integer"); }
	}

	try { validateRequest(request); }
	catch (...) { throw; }

	return request;
}

SimulationRequest newSimulationRequest(std::string npop, std::string prevalence,
	std::string nsnp, std::string me, std::string sde, std::string oddsRatio,
	std::string lowerMAF, std::string upperMAF, std::vector<SimulationRequestGroup> groups,
	bool useCommonTest, bool useBootstrap, std::string nboot) {

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
	try { request.upperMAF = std::abs(std::stod(upperMAF)); }
	catch (...) { throw std::invalid_argument("MAF upper bound,decimal"); }
	try { request.lowerMAF = std::abs(std::stod(lowerMAF)); }
	catch (...) { throw std::invalid_argument("MAF lower bound,decimal"); }

	request.useCommonTest = useCommonTest;
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

std::vector<std::vector<double>> startVikNGS(Request req) {

	VectorXd Y, G; MatrixXd X, Z, P;
	std::map<int, int> readGroup;
	std::vector<std::vector<int>> interval;

	bool valid = parseAndFilter(req.vcfDir, req.sampleDir, req.bedDir,
		req.highLowCutOff, req.collapseCoding, req.collapseExon,
		req.missingThreshold, req.onlySNPs, req.mustPASS, req.mafCutoff, req.useCommonTest,
		X, Y, Z, G, readGroup, P, interval);

	bool useCovariates = Z.rows() > 0;

	std::vector<std::vector<double>> pvals;


	if (!valid) {
		//TODO 
		return pvals;
	}

	if (req.useCommonTest) {

		std::vector<double> pval;

		if (req.useBootstrap && useCovariates)
			pval = runCommonTest(X, Y, Z, G, readGroup, P, req.nboot);
		else if(req.useBootstrap && !useCovariates)
			pval = runCommonTest(X, Y, G, readGroup, P, req.nboot);
		else if (!req.useBootstrap && useCovariates)
			pval = runCommonTest(X, Y, Z, G, readGroup, P);
		else
			pval = runCommonTest(X, Y, G, readGroup, P);

		for (int i = 0; i < pval.size(); i++) {
			std::vector<double> p;
			p.push_back(pval[i]);
			pvals.push_back(p);
		}
	}
	else {
		if (useCovariates)
			pvals = runRareTest(X, Y, Z, G, readGroup, P, req.nboot);
		else
			pvals = runRareTest(X, Y, G, readGroup, P, req.nboot);
	}

	std::string outputDir = "C:/Users/Scott/Desktop/out.txt";
	outputPvals(pvals, outputDir);

	return pvals;
}

std::vector<std::vector<double>> startSimulation(SimulationRequest req) {

	VectorXd Y, G; MatrixXd X, P;
	std::map<int, int> readGroup;

	simulate(req, X, Y, G, readGroup, P);

	std::vector<std::vector<double>> pvals;

	if (req.useCommonTest){
		
		std::vector<double> pval;

		if(req.useBootstrap)
			pval = runCommonTest(X, Y, G, readGroup, P, req.nboot);
		else
			pval = runCommonTest(X, Y, G, readGroup, P);

		for (int i = 0; i < pval.size(); i++) {
			std::vector<double> p;
			p.push_back(pval[i]);
			pvals.push_back(p);
		}
	}
	else
		pvals = runRareTest(X, Y, G, readGroup, P, req.nboot);

    std::string outputDir = "C:/Users/Scott/Desktop/out.txt";
    outputPvals(pvals, outputDir);

    return pvals;
}

SimulationRequest testSimulationRequest() {
	std::vector<SimulationRequestGroup> groups;
	groups.push_back(newSimulationRequestGroup(0, "300", "control", "low", "10", "5"));
	groups.push_back(newSimulationRequestGroup(0, "200", "case", "high", "40", "7"));
	groups.push_back(newSimulationRequestGroup(0, "100", "control", "high", "50", "6"));

	return newSimulationRequest("3000", "0.1", "100", "0.01", "0.025", "1.4", "0.1", "0.5", groups, true, false, 0);
}

/*
int main() {

    //TODO: take as input from command line
    //---------------------------------------
    bool simulation = true;
    bool common = true;

    int highLowCutOff = 30;
    bool collapseCoding = false;
    bool collapseExon = true;

    //filtering paramaters
    double mafCutoff = 0.05;
    double missingThreshold = 0.2;
    bool onlySNPs = true;
    bool mustPASS = true;

    bool rvs = true;

    ///input files
    std::string vcfDir = "C:/Users/Scott/Desktop/RVS-master/example/example_1000snps.vcf";
    std::string infoDir = "C:/Users/Scott/Desktop/RVS-master/example/sampleInfo.txt";
    //std::string bedDir = "";
    std::string bedDir = "C:/Users/Scott/Desktop/RVS-master/example/chr11.bed";

    std::string outputDir = "C:/Users/Scott/Desktop/out.txt";

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
        bool valid = parseAndFilter(vcfDir, infoDir, bedDir,
            highLowCutOff, collapseCoding, collapseExon,
            missingThreshold, onlySNPs, mustPASS, mafCutoff, common,
            X, Y, Z, G, readGroup, P, interval);
        if (!valid) {
            return 0;
        }


        //generateForR(X, Y, Z, G, P, readGroup);

        if (common) {
            //std::vector<double> pvals = runCommonTest(X, Y, G, readGroup, P, 1000);
            std::vector<double> pvals = runCommonTest(X, Y, Z, G, readGroup, P);
            //std::vector<double> pvals = runCommonTest(X, Y, Z, G, readGroup, P, 1000, true);


            std::cout << "Common Test p-values\n";
            for (size_t i = 0; i < pvals.size(); i++) {
                std::cout << pvals[i];
                std::cout << '\n';
            }
        }
        else {

            std::vector<std::vector<double>> pval = runRareTest(X, Y, G, readGroup, P, 20000, true);
            //std::vector<std::vector<double>> pval = runRareTest(X, Y, G, readGroup, P, Z, 20000, true);

            std::cout << "Rare Test p-values\n";
            std::cout << pval[0][0];
            std::cout << '\t';
            std::cout << pval[0][1];

        }
    }


    //keep console open while debugging
    //TODO: be sure to remove eventually!
    std::cout << "\ndone...>";

    while (true) {}
    return 0;
}

*/
