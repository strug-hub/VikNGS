#pragma once
#include "../vikNGS.h"
#include "../Math/Math.h"

class TestObject {

    MatrixXd X;
    VectorXd Y;
    MatrixXd Z;
    VectorXi G;
    std::map<int, Depth>& d;
    MatrixXd P;
    Family family;

    VectorXd Ycenter;
    VectorXd Ycenter_original;

    VectorXd MU;

    MatrixXd Xboot;
    VectorXd Yboot;
    MatrixXd Zboot;

    inline void calculateYCenter() {

        if(hasCovariates()){
            VectorXd beta = getBeta(Y, Z, family);
            this->MU = fitModel(beta, Z, family);
        }
        else{
            double ybar = Y.mean();
            this->MU = VectorXd::Constant(this->Y.rows(), ybar);
        }

        Ycenter = Y - MU;
    }

public:

    TestObject(MatrixXd& genotypes, VectorXd& phenotypes, MatrixXd& covariates,
               MatrixXd& frequency, Family distribution, VectorXi& groups, std::map<int, Depth>& depths, bool rareVariant) :
        X(genotypes), Y(phenotypes), Z(covariates), G(groups), d(depths), P(frequency), family(distribution) {




        //Filter NAN
        VectorXi toRemove = whereNAN(Y);
        if(!rareVariant)
            toRemove = toRemove + whereNAN(X);
        if(hasCovariates())
            toRemove = toRemove + whereNAN(Z);

        X = extractRows(X, toRemove, 0);
        Y = extractRows(Y, toRemove, 0);
        if(hasCovariates())
            Z = extractRows(Z, toRemove, 0);
        G = extractRows(G, toRemove, 0);

        //Replace NAN with 0 if rare
        if(rareVariant)
            X = replaceNAN(X, 0);

        calculateYCenter();

        bootstrapped = false;
        groupVectorCache = false;
        XcenterCache = false;
    }

    inline bool hasCovariates() { return this->Z.cols() > 1; }
    inline double getRobustVar(int i=0){
        return calcRobustVar(P(i,1), P(i,2));
    }
    inline VectorXd robustVarVector(){
        VectorXd robustVar(P.rows());
        for (int i = 0; i < P.rows(); i++)
            robustVar[i] = sqrt(getRobustVar(i));

        return robustVar;
    }
    inline VectorXd mafWeightVector(){
        VectorXd maf(P.rows());
        for (int i = 0; i < P.rows(); i++){
            double m = P(i,0) + 0.5*P(i,1);
            maf[i] = 1/sqrt(m * (1-m));
        }

        return maf;
    }

    inline MatrixXd* getX(){ return (bootstrapped) ? &Xboot : &X; }
    inline VectorXd getX(int i){ return (bootstrapped) ? Xboot.col(i) : X.col(i); }
    inline VectorXd* getY(){ return (bootstrapped) ? &Yboot : &Y; }
    inline MatrixXd* getZ(){ return (bootstrapped) ? &Zboot : &Z; }
    inline VectorXi* getG(){ return &G; }
    inline VectorXd* getMU(){ return &MU; }
    inline VectorXd* getYcenter(){ return &Ycenter; }

    inline std::map<int, Depth>* getDepths(){ return &d; }

    inline void bootstrap(Test& test, Family family) {
        calculateXcenter();

       if(test.isExpectedGenotypes()){
           calculateGroupVector();

           if(family == Family::NORMAL)
               normalBootstrap();
           if(family == Family::BINOMIAL)
               binomialBootstrap();
       }
       else
           permute();

       bootstrapped = true;
       calculateYCenterBoot();
   }


//bootstrap functions
private:

    bool bootstrapped;

    inline void permute() {
        if(!bootstrapped)
           Yboot = Y;
        Xboot = shuffleColumnwiseWithoutReplacement(Xcenter);
        if(hasCovariates())
           Zboot = shuffleColumnwiseWithoutReplacement(Z);
    }

    void binomialBootstrap() {
        if(!bootstrapped)
            Yboot = Y;
        Xboot = groupwiseShuffleWithReplacement(Xcenter, G, groupVector);

        if(hasCovariates())
            Zboot = groupwiseShuffleWithReplacement(Z, G, groupVector);
    }

    void normalBootstrap() {

       if(hasCovariates()){

           if(!bootstrapped)
               Xboot = X;
               Zboot = Z;
               Ycenter_original = Ycenter;

           VectorXd residuals = groupwiseShuffleWithoutReplacement(Ycenter_original, G, groupVector);
           Yboot = Y - Ycenter_original + residuals;

       }
       else{
           Yboot = groupwiseShuffleWithoutReplacement(Y, G, groupVector);
           if(!bootstrapped)
               Xboot = X;
       }
    }

    bool groupVectorCache;
    std::map<int, std::vector<int>> groupVector;

    inline void calculateGroupVector() {

        if(groupVectorCache)
            return;

        for(int i = 0; i < G.rows(); i++){

            if(groupVector.count(G[i]) < 1){
                std::vector<int> v;
                v.push_back(i);
                groupVector[G[i]] = v;
            }
            else
                groupVector[G[i]].push_back(i);
        }

        groupVectorCache = true;
    }

    bool XcenterCache;
    MatrixXd Xcenter;
    inline void calculateXcenter() {
        if(XcenterCache)
            return;

        Xcenter = subtractGroupMean(X, G);
    }

    inline void calculateYCenterBoot() {

        if(hasCovariates()){
            VectorXd beta = getBeta(Yboot, Zboot, family);
            this->MU = fitModel(beta, Zboot, family);
        }
        else{
            double ybar = Yboot.mean();
            this->MU = VectorXd::Constant(Yboot.rows(), ybar);
        }

        Ycenter = Yboot - MU;
    }
};


   /*

    std::vector<MatrixXd> h;
    std::vector<VectorXd> hdiag;
    VectorXd Yvar;
    std::vector<VectorXd> yvar;

   void constructHatMatrix(Family family, VectorXd &X, VectorXd &Y, MatrixXd &Z, VectorXd &G, std::vector<VectorXd> &yhat){

       if(family == Family::NORMAL){

           H = calculateHatMatrix(Z);
           double constant = variance(Y);

           for (size_t i = 0; i < x.size(); i++)
                   this->yvar.push_back(MatrixXd::Constant(x[i].rows(), 1, constant));
       }
       else if(family == Family::BINOMIAL){
           VectorXd MU_ = MU.array() * (VectorXd::Ones(MU.rows()) - MU).array();
           VectorXd MU_sqrt = mu_.array().sqrt();

           MatrixXd W = mu_.asDiagonal();
           MatrixXd sqrtW = mu_sqrt.asDiagonal();

           H = calculateHatMatrix(Z_filtered, W, sqrtW);

           for (int i = 0; i < size(); i++)
               this->yvar.push_back(extractRows(mu_, G_filtered, i).transpose());
       }

       VectorXd Hdiag = H.diagonal();

       for (int i = 0; i < size(); i++){
           this->h.push_back(extractRows(H, G_filtered, i).transpose());
           this->hdiag.push_back(extractRows(Hdiag, G_filtered, i));
       }

       this->Yvar = concatenate(this->yvar);
   }

   */
