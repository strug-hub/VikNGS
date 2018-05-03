#include "../Math/MathHelper.h"
#include "../Request.h"
#include "../TestInput.h"
#include "RareTestObject.h"
#include <future>

double getTestStatistic(std::vector<RareTestObject> &t, bool rvs, std::string test) {
	int nsnp = t.size();
	int i, j;

	MatrixXd diagS = MatrixXd::Constant(nsnp, nsnp, 0);

	VectorXd var(nsnp);
	VectorXd ym_hrd(nsnp);
	VectorXd ym_lrd(nsnp);

	for (i = 0; i < nsnp; i++) {
		var[i] = sqrt(t[i].getRobustVar());
		ym_hrd[i] = t[i].getYm(1);
		ym_lrd[i] = t[i].getYm(0);
	}

	MatrixXd diagVar = var.asDiagonal();
	MatrixXd diagYm_hrd = ym_hrd.asDiagonal();
	MatrixXd diagYm_lrd = ym_lrd.asDiagonal();

	for (i = 0; i < t[0].size(); i++) {
		VectorXd x_0 = t[0].getX(i);
		MatrixXd x(x_0.rows(), nsnp);
		x.col(0) = x_0;

		for (j = 1; j < nsnp; j++)
			x.col(j) = t[j].getX(i);

		if (t[0].isHRG(i)) {
			if (rvs) {
				MatrixXd var_hrd = diagVar.transpose() * correlation(x) * diagVar;
				diagS += diagYm_hrd * var_hrd * diagYm_hrd;
			}
			else
				diagS += diagYm_hrd * covariance(x) * diagYm_hrd;
		}
		else
			diagS += diagYm_lrd * covariance(x) * diagYm_lrd;
	}

	VectorXd score(nsnp);
	for (i = 0; i < nsnp; i++)
		score[i] = t[i].getScore();

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

double rareTest(std::vector<RareTestObject> &t, int nboot, bool stopEarly, bool rvs, std::string test) {
	int i, j, h;
	int nsnp = t.size();

	double tobs = getTestStatistic(t, rvs, test);

	//start bootstrapping!

	double bootCount = 0;
	double tcount = 0;
	double tsamp;

	for (h = 0; h < nboot; h++) {
		for (i = 0; i < nsnp; i++)
			t[i].bootstrap();

		tsamp = getTestStatistic(t, false, test);

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

	return{ (tcount + 1) / (bootCount + 1)};
}

std::vector<double> runRareTestParallel(Request &req, TestInput &input, int threadID, int nthreads) {
  
  MatrixXd X = input.X;
  VectorXd Y = input.Y;
  MatrixXd Z = input.Z;
  
  VectorXd G = input.G;
  std::map<int, int> readGroup = input.readGroup;
  MatrixXd P = input.P;
  
  //todo calpha cannot go in parallel???
  std::string test = req.rareTest;

  int i, j, k;
  std::vector<double> pvals;
  
  std::vector<MatrixXd> x;
  std::vector<VectorXd> y;
  std::vector<MatrixXd> z;
  std::vector<int> rd;
  
  int ngroups = 1 + (int)G.maxCoeff();
  bool hasCovariates = input.hasCovariates();
  
  for (i = 0; i < ngroups; i++) {
    x.push_back(extractRows(X, G, i));
    y.push_back(extractRows(Y, G, i));
    
    if (hasCovariates)
      z.push_back(extractRows(Z, G, i));
    
    rd.push_back(readGroup[i]);
  }

  for (k = threadID; k < input.interval.size(); k+=nthreads) {

    if (k % (25) == 0) {
      std::cout << std::endl << k;
      std::cout << "/" << input.interval.size();
      std::cout << " calculated." << std::endl;

    }
    
    std::vector<RareTestObject> t;

      
    std::cout << "Collasping: ";
    for (i = 0; i < input.interval[k].size(); i++) {
    	std::cout << input.interval[k][i] << " ";
    }
    std::cout << "\n";

    try{
        for (i = 0; i < input.interval[k].size(); i++) {

          int index = input.interval[k][i];

          std::vector<VectorXd> x_i;
          for (j = 0; j < ngroups; j++)
            x_i.push_back(x[j].col(index));

          VectorXd X_i = X.col(index);
          VectorXd P_i = P.row(index);

          if (hasCovariates) {
            RareTestObject t_i(X_i, Y, Z, x_i, y, z, rd, P_i);
            t.push_back(t_i);
          }
          else{
            RareTestObject t_i(x_i, y, rd, P_i);
            t.push_back(t_i);
          }
        }

        pvals.push_back(rareTest(t, req.nboot, req.stopEarly, req.rvs, test));

    }catch(variant_exception){

            //todo: printinfo for collapse
            printWarning("Problem occured while running test on a variant block (...). Assigning a p-value of 1");
            pvals.push_back(1);
        }
  }
  
  return pvals;
}


std::vector<Variant> runRareTest(Request req, TestInput input) {
  int nthreads = req.nthreads;
  std::vector<std::future<std::vector<double>>> pvalsFuture;
  
  if (nthreads <= 1){
    std::vector<double> pvals = runRareTestParallel(req, input, 0, 1);
    for (int i = 0; i < pvals.size(); i++)
        input.variants[i].pvalue = pvals[i];
    return input.variants;
  }
  
  for (int i = 0; i < nthreads; i++) {
    pvalsFuture.push_back(std::async(std::launch::async,
                                     [&,i] { return runRareTestParallel(req, input, i, nthreads); } ));
  }
  
  std::vector<std::vector<double>> pvalsFromThread;
  
  for (int i = 0; i < nthreads; i++) 
    pvalsFromThread.push_back(pvalsFuture[i].get());
  
  std::vector<double> pvals;
  
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
