#pragma once
#include "vikNGS.h"
#include "Math/Math.h"
#include "Eigen/Dense"
#include "Interval.h"

#include <string>
#include <vector>
#include <map>

using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::MatrixXd;

struct Variant {
private:
    std::map<Genotype, VectorXd> genotypes;
    std::map<Genotype, Vector3d> P;

    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::string alt;

    Filter filter;

public:

    Variant(std::string chromosome, int position, std::string unique_id, std::string reference, std::string alternative) :
        chrom(chromosome), pos(position), id(unique_id), ref(reference), alt(alternative) {
        filter = Filter::VALID;
    }

    Variant() {
        filter = Filter::INVALID;
    }
    ~Variant() { }

    inline void setExpectedGenotypes(std::vector<Vector3d> &likelihoods) {
        P[Genotype::EXPECTED] = calculateGenotypeFrequencies(likelihoods);
        this->genotypes[Genotype::EXPECTED] = calculateExpectedGenotypes(likelihoods, P[Genotype::EXPECTED]);
    }
    inline void setCallGenotypes(std::vector<Vector3d> &likelihoods) {
        P[Genotype::CALL] = calculateGenotypeFrequencies(likelihoods);
        this->genotypes[Genotype::CALL] = calculateGenotypeCalls(likelihoods, P[Genotype::CALL]);
    }
    inline void setTrueGenotypes(VectorXd &gt) {
        P[Genotype::TRUE] = calculateGenotypeFrequencies(gt);
        this->genotypes[Genotype::TRUE] = gt;
    }
    inline void setVCFCallGenotypes(VectorXd &gt) {
        P[Genotype::VCF_CALL] = calculateGenotypeFrequencies(gt);
        this->genotypes[Genotype::VCF_CALL] = gt;
    }

    inline void setFilter(Filter f) { this->filter = f; }
    inline std::string getChromosome() { return this->chrom; }
    inline int getPosition() { return this->pos; }

    inline std::vector<Genotype> getAllGenotypes(){
        std::vector<Genotype> all;
        for(auto const& gt : this->genotypes)
            all.push_back(gt.first);
        return all;
    }
    inline VectorXd getGenotype(Genotype gt) { return genotypes[gt]; }
    inline Vector3d getP(Genotype gt) { return P[gt]; }

    inline bool isValid() { return this->filter == Filter::VALID; }

    inline std::string toString() {
        std::string t = "\t";
        return chrom + t + std::to_string(pos) + t + ref + t + alt;
    }
    inline std::string getInfo() {
        std::string t = "\t";
        return "chr:" + chrom + " pos:" + std::to_string(pos) + " ref:" + ref + " alt:" + alt;
    }
    inline bool operator < (Variant& v) {
        if (this->chrom == v.getChromosome())
            return this->pos < v.getPosition();
        else
            return this->chrom < v.getChromosome();
    }

};

inline bool variantCompare(Variant lhs, Variant rhs) { return lhs < rhs; }

struct VariantSet{
private:
    std::vector<Variant> variants;
    std::map<Test, double> pval;
    int nvalid = 0;
    Interval *interval;
public:
    VariantSet(Variant & v) { variants.push_back(v); }
    VariantSet() { }

    inline void setInterval(Interval * inv) { interval = inv; }
    inline bool isIn(Variant &variant) { return interval->isIn(variant.getChromosome(), variant.getPosition()); }

    inline void addVariant(Variant &variant) {
        variants.push_back(variant);
        if(variant.isValid()) nvalid++;
    }
    inline void addPval(Test test, double p) { this->pval[test] = p; }

    inline double getPval(Test test) { return this->pval[test]; }
    inline int nPvals() { return this->pval.size(); }
    inline int size(){ return variants.size(); }
    inline int validSize(){ return variants.size(); }

};









    double trueMaf;


    std::vector<int> readDepths;
    std::vector<std::vector<int>> baseCalls;
    std::vector<std::string> vcfCalls;
    std::string format;
    std::vector<std::string> columnUsed;









    inline void reduceSize() {
        likelihood.clear();
        VectorXd empty;
        P = empty;
        expectedGenotype = empty;
    }





    inline MatrixXd getLikelihoodMatrix() {

        MatrixXd M(likelihood.size(), 3);
        for(int i = 0; i < likelihood.size(); i++){
            M(i,0)=likelihood[i].L00;
            M(i,1)=likelihood[i].L01;
            M(i,2)=likelihood[i].L11;
        }

        return M;
    }

    //Uses EM algorithm to estimate the genotype frequencies in the sample.
    inline void calculateGenotypeFrequency() {

        MatrixXd M = getLikelihoodMatrix();

        double p = 0.15;
        double q = 0.15;
        double qn = 1;
        double pn = 0;
        double dn = 0;
        double d = 0;

        int k = 0;

        while (pow((pn - p), 2) + pow((qn - q), 2) > 0.000001) {

            d = 1 - p - q;
            Vector3d v = { p, q, d };
            VectorXd pD = M * v;
            VectorXd Ep = M.col(0).array() * (p / pD.array()).array();
            VectorXd Eq = M.col(1).array() * (q / pD.array()).array();
            pn = p;
            qn = q;
            dn = 1 - q - p;
            p = Ep.sum() / Ep.rows();
            q = Eq.sum() / Eq.rows();

            k++;
            if (k == 1000)
                break;
        }

        VectorXd freq(3);
        freq[0] = std::max(0.0, p);
        freq[1] = std::max(0.0, q);
        freq[2] = std::max(0.0, 1 - p - q);

        P=freq;
    }

};



