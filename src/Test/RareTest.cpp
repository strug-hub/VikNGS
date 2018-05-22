#include "../Math/MathHelper.h"
#include "../Request.h"
#include "../TestInput.h"
#include "RareTestObject.h"
#include "RareTestCollapseObject.h"

#include <future>

double getTestStatistic(RareTestCollapseObject &collapse, bool rvs, std::string test) {

    MatrixXd diagS = collapse.getVariance(rvs);
    VectorXd score = collapse.getScore();

	VectorXd e = diagS.eigenvalues().real();
	std::vector<double> eigenvals(e.data(), e.data() + e.size());

	double tstat = 0;

	//CAST -> linear test
	if(test == "cast")
		tstat = pnorm(score.sum() / sqrt(diagS.sum()));
	//C-alpha -> quadratic test
	else if(test == "calpha")
		tstat = qfc(eigenvals, score.array().pow(2).sum(), score.rows());

	return tstat;
}

double rareTest(RareTestCollapseObject &collapse, int nboot, bool stopEarly, bool rvs, std::string test) {
    int h;
    double tobs = getTestStatistic(collapse, rvs, test);

	//start bootstrapping!

	double bootCount = 0;
	double tcount = 0;
	double tsamp;

	for (h = 0; h < nboot; h++) {
        collapse.bootstrap();

        tsamp = getTestStatistic(collapse, false, test);

		if (tsamp <= tobs)
			tcount++;

		bootCount++;

		if (stopEarly && bootCount > 100) {
			double pstar = 5 / ((bootCount + 0.0737224)*(1 + 0.0737224));
			if (tcount / h > pstar) {
				std::cout << "early stop at ";
				std::cout << bootCount;
				std::cout << "\n ";
				break;
			}
		}
	}

    //todo: pvalue is 0???
    return (tcount + 1) / (bootCount + 1) ;
}

std::vector<double> runRareTestParallel(std::vector<MatrixXd> x, std::vector<VectorXd> y, std::vector<MatrixXd> z,
                                        std::vector<int> rd, std::vector<VectorXd> ycenter,
                                        Request &req, TestInput &input, int threadID, int nthreads) {
  
    MatrixXd P = input.P;
    int ngroups = 1 + (int)input.G.maxCoeff();

    //todo calpha cannot go in parallel???
    std::string test = req.rareTest;
    std::vector<double> pvals;

  for (int k = threadID; k < input.interval.size(); k+=nthreads) {

    if (k % (25) == 0) {
      std::cout << std::endl << k;
      std::cout << "/" << input.interval.size();
      std::cout << " calculated." << std::endl;

    }
    
    RareTestCollapseObject collapse(y, z, ycenter, input.hasCovariates(), input.family);

    std::cout << "Collasping: ";
    for (int i = 0; i < input.interval[k].size(); i++) {
    	std::cout << input.interval[k][i] << " ";
    }
    std::cout << "\n";

    try{
        for (int i = 0; i < input.interval[k].size(); i++) {

          int index = input.interval[k][i];

          std::vector<VectorXd> x_i;
          for (int j = 0; j < ngroups; j++)
            x_i.push_back(x[j].col(index));

          //VectorXd X_i = X.col(index);
          VectorXd p = P.row(index);
          collapse.addVariant(x_i, rd, p);

        }

        pvals.push_back(rareTest(collapse, req.nboot, req.stopEarly, req.rvs, test));

    }catch(variant_exception){

            //todo: printinfo for collapse
            printWarning("Problem occured while running test on a variant block (...). Assigning a p-value of 1");
            pvals.push_back(1);
        }
  }
  
  return pvals;
}


std::vector<Variant> runRareTest(Request req, TestInput input) {

    //preprocessing
    //------------------

    VectorXd Y = input.Y;
    MatrixXd Z = input.Z;
    bool covariates = input.hasCovariates();
    MatrixXd X = input.X;
    VectorXd G = input.G;
    VectorXd toRemove;

    if(covariates){
        toRemove = whereNAN(Y, Z);
        Z = extractRows(Z, toRemove, 0);
    }
    else
        toRemove = whereNAN(Y);

    Y = extractRows(Y, toRemove, 0);
    G = extractRows(G, toRemove, 0);
    X = extractRows(X, toRemove, 0);
    X = replaceNAN(X,  0);

    std::map<int, int> readGroup = input.readGroup;

    std::vector<MatrixXd> x;
    std::vector<VectorXd> y;
    std::vector<MatrixXd> z;
    std::vector<int> rd;

    int ngroups = 1 + (int)G.maxCoeff();

    for (int i = 0; i < ngroups; i++) {
      x.push_back(extractRows(X, G, i));
      y.push_back(extractRows(Y, G, i));

      if (covariates)
        z.push_back(extractRows(Z, G, i));

      rd.push_back(readGroup[i]);
    }

    std::vector<VectorXd> ycenter;

  try{
      if(covariates){
          VectorXd beta = getBeta(Y, Z, input.family);
          ycenter = fitModel(beta, y, z, input.family);
      }
      else{
          double ybar = average(y);
          for (int i = 0; i < y.size(); i++)
              ycenter.push_back(y[i].array() - ybar);
      }

  }catch(...){
      //todo when will this error be caught? correlated response variables?
      printWarning("Error occured while regressing covariates on response variable.");
  }

  //------------------

  int nthreads = req.nthreads;
  std::vector<std::future<std::vector<double>>> pvalsFuture;
  
  if (nthreads <= 1){
    std::vector<double> pvals = runRareTestParallel(x, y, z, rd, ycenter, req, input, 0, 1);
    for (int i = 0; i < pvals.size(); i++)
        input.variants[i].pvalue = pvals[i];
    return input.variants;
  }
  
  for (int i = 0; i < nthreads; i++) {
    pvalsFuture.push_back(std::async(std::launch::async,
                                     [&,i] { return runRareTestParallel(x, y, z, rd, ycenter, req, input, i, nthreads); } ));
  }
  
  std::vector<std::vector<double>> pvalsFromThread;
  
  for (int i = 0; i < nthreads; i++) 
    pvalsFromThread.push_back(pvalsFuture[i].get());
    
  int i = 0;
  int index = 0;
  while (true) {
    
    bool flag = false;
    for (int j = 0; j < nthreads; j++)
      if (pvalsFromThread[j].size() > i) {
        input.variants[index].pvalue = pvalsFromThread[j][i];
        index++;
        flag = true;
      }
      
      i++;
      
      if (!flag)
        break;
  }
  
  return input.variants;
}
