#pragma once
#include "Enum/Family.h"
#include "Enum/Depth.h"

#include <map>
#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

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
            if( !(std::abs(Y[i]) < epsilon || std::abs(Y[i] - 1) < epsilon)){
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
    inline void setFamily(Family fam){ family = fam; }

    inline bool hasCovariates() { return Z.rows() > 0 && Z.cols() > 0; }
    inline int ncov() { return Z.cols() -1; }
    inline int nsamp() { return Y.rows(); }
    inline int ngroup() { return 1 + G.maxCoeff(); }
    inline std::map<int, Depth> getGroupDepthMap(){ return groupDepth; }

};

