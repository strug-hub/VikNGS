#pragma once
#include "vikNGS.h"
#include <vector>
#include <map>

#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

struct SampleInfo {
private:
    VectorXd Y;
    VectorXi G;
    MatrixXd Z;
    std::map<int, Depth> groupDepth;
    Family family;

    void determineFamily() {
        double epsilon = 1e-8;
        //if a value not 0 or 1 is found, assume quantitative data
        for(int i = 0; i < Y.rows(); i++){
            if( std::abs(Y[i]) > epsilon || std::abs(Y[i] - 1) > epsilon){
                family=Family::NORMAL;
                return;
            }
        }

        family=Family::BINOMIAL;
    }

public:

    inline void setY(VectorXd Y){ this->Y = Y; determineFamily(); }
    inline void setZ(MatrixXd Z){ this->Z = Z; }
    inline void setG(VectorXi G){ this->G = G; }
    inline void setGroupDepthMap(std::map<int, Depth> groupDepth){ this->groupDepth = groupDepth; }

    inline VectorXi getG(){ return G; }
    inline VectorXd getY(){ return Y; }
    inline MatrixXd getZ(){ return Z; }
    inline Family getFamily(){ return family; }

    inline bool hasCovariates() { return Z.rows() > 0 && Z.cols() > 0; }
    inline int ncov() { return Z.cols() -1; }
    inline int nsamp() { return Y.rows(); }
    inline int ngroup() { return 1 + G.maxCoeff(); }

};




/*


        inline MatrixXd getP(){
            MatrixXd P(variants.size(), 3);
            for (int i = 0; i < variants.size(); i++)
                    P.row(i) = variants[i].P;
            return P;
        }

        inline MatrixXd getX(bool useRegular, bool useTrueGenotypes){

            MatrixXd X(variants[0].likelihood.size(), variants.size());

            if(!useRegular){
            for (int i = 0; i < variants.size(); i++)
                    X.col(i) = variants[i].expectedGenotype;
            }
            else if(useRegular && useTrueGenotypes){
                for (int i = 0; i < variants.size(); i++)
                    X.col(i) = variants[i].trueGenotype;
            }
            else{
                for (int i = 0; i < variants.size(); i++)
                    X.col(i) = variants[i].genotypeCalls;
            }

            return X;
        }

        inline void addCollapse(std::vector<std::vector<int>> c) {
                collapse = c;
        }
};


//Creates a TestInput object.

inline SampleInfo buildTestInput(VectorXd &Y, MatrixXd &Z, VectorXd &G, std::map<int, int> &readGroup,
                                std::vector<Interval> intervals, std::string family) {

        SampleInfo input;
	input.Y = Y;
	input.Z = Z;
	input.G = G;
	input.readGroup = readGroup;
        input.family = family;
        input.intervals = intervals;
	return input;
}


//Creates a TestInput object.

inline SampleInfo buildTestInput(VectorXd &Y, MatrixXd &Z, VectorXd &G,
        std::map<int, int> &readGroup, std::vector<Variant> &variants, std::string family) {

        SampleInfo input;
        input.Y = Y;
        input.Z = Z;
        input.G = G;
        input.readGroup = readGroup;
        input.variants = variants;
        input.family = family;
        return input;
}

inline SampleInfo addVariants(SampleInfo t, std::vector<Variant> &variants, std::vector<std::vector<int>> &collapse) {

        t.variants = variants;
        t.collapse = collapse;
        return t;
}


*/
