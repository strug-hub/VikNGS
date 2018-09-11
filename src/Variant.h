#pragma once
#include "Enums.h"
#include "./Math/Math.h"
#include "Interval.h"
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;

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
    double trueMaf;

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
    inline void setTrueGenotypes(VectorXd& gt) {
        P[Genotype::TRUE] = calculateGenotypeFrequencies(gt);
        this->genotypes[Genotype::TRUE] = gt;
    }
    inline void setVCFCallGenotypes(VectorXd& gt) {
        P[Genotype::VCF_CALL] = calculateGenotypeFrequencies(gt);
        this->genotypes[Genotype::VCF_CALL] = gt;
    }

    inline void setFilter(Filter f) { this->filter = f; }
    inline void setTrueMaf(double maf) { this->trueMaf = maf; }

    inline std::string getChromosome() { return this->chrom; }
    inline std::string getRef() { return this->ref; }
    inline std::string getAlt() { return this->alt; }

    inline int getPosition() { return this->pos; }

    inline std::vector<Genotype> getAllGenotypes(){
        std::vector<Genotype> all;
        for(auto const& gt : this->genotypes)
            all.push_back(gt.first);
        return all;
    }
    //inline VectorXd getGenotype(Genotype gt) { return genotypes[gt]; }
    inline VectorXd* getGenotype(Genotype gt) { return &genotypes[gt]; }

    inline Vector3d* getP(Genotype gt) { return &P[gt]; }

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

    inline void reduceSize() {
        genotypes.clear();
    }


};

inline bool variantCompare(Variant lhs, Variant rhs) { return lhs < rhs; }

struct VariantSet{
private:
    std::vector<Variant> variants;
    std::vector<double> pval;
    int nvalid = 0;
    Interval *interval;
public:
    VariantSet(Variant& v) { variants.push_back(v); }
    VariantSet() { }

    inline void setInterval(Interval * inv) { interval = inv; }
    inline bool isIn(Variant &variant) { return interval->isIn(variant.getChromosome(), variant.getPosition()); }

    inline void addVariant(Variant &variant) {
        variants.push_back(variant);
        if(variant.isValid()) nvalid++;
    }

    inline std::vector<Variant>* getVariants() { return &variants; }
    inline std::string getChromosome() { return (variants.size() < 1) ? "na" : variants[0].getChromosome(); }
    inline int getMinPos() { return (variants.size() < 1) ? -1 : variants[0].getPosition(); }
    inline int getMaxPos() { return (variants.size() < 1) ? -1 : variants.back().getPosition(); }
    inline double getMidPos() { return (variants.size() < 1) ? -1 :
                                                            (variants[0].getPosition() + variants.back().getPosition())/2.0; }

    inline void addPval(double p) { pval.push_back(p); }

    inline double getPval(size_t i) { return pval[i]; }
    inline int nPvals() { return static_cast<int>(this->pval.size()); }
    inline int size(){ return static_cast<int>(variants.size()); }
    inline int validSize(){ return nvalid; }
    inline bool isValid(){ return nvalid > 0; }

    inline MatrixXd getX(Genotype gt){
        std::vector<VectorXd*> x;
        for (size_t i = 0; i < static_cast<size_t>(validSize()); i++)
            if(variants[i].isValid())
                x.push_back(variants[i].getGenotype(gt));

        if(x.size() < 1)
            return MatrixXd();

        MatrixXd X(x[0]->rows(), x.size());
        for (size_t i = 0; i < x.size(); i++)
            X.col(static_cast<int>(i)) = *x[i];

        return X;
    }

    inline MatrixXd getP(Genotype gt){
        std::vector<Vector3d*> p;
        for (size_t i = 0; i < static_cast<size_t>(validSize()); i++)
            if(variants[i].isValid())
                p.push_back(variants[i].getP(gt));

        if(p.size() < 1)
            return MatrixXd();

        MatrixXd P(p.size(), 3);
        for (size_t i = 0; i < p.size(); i++)
            P.row(static_cast<int>(i)) = *p[i];

        return P;
    }
};
