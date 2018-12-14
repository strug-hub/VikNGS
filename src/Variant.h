#pragma once
#include "./Math/Math.h"
#include "Interval.h"

#include <string>
#include <iostream>


struct Variant {
private:
    std::map<GenotypeSource, VectorXd> genotypes;
    std::map<GenotypeSource, Vector3d> P;

    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::string alt;

    Filter filter;
    bool isShrunk;

public:

    Variant(std::string chromosome, int position, std::string unique_id, std::string reference, std::string alternative) :
        chrom(chromosome), pos(position), id(unique_id), ref(reference), alt(alternative) {
        filter = Filter::VALID;
        isShrunk = false;
    }

    Variant() {
        filter = Filter::INVALID;
        isShrunk = false;
    }
    ~Variant() { }

    inline void setExpectedGenotypes(std::vector<Vector3d> &likelihoods) {
        P[GenotypeSource::EXPECTED] = calculateGenotypeFrequencies(likelihoods);
        this->genotypes[GenotypeSource::EXPECTED] = calculateExpectedGenotypes(likelihoods, P[GenotypeSource::EXPECTED]);
    }
    inline void setCallGenotypes(std::vector<Vector3d> &likelihoods) {
        P[GenotypeSource::CALL] = calculateGenotypeFrequencies(likelihoods);
        this->genotypes[GenotypeSource::CALL] = calculateGenotypeCalls(likelihoods, P[GenotypeSource::CALL]);
    }
    inline void setTrueGenotypes(VectorXd& gt) {
        P[GenotypeSource::TRUE] = calculateGenotypeFrequencies(gt);
        this->genotypes[GenotypeSource::TRUE] = gt;
    }
    inline void setVCFCallGenotypes(VectorXd& gt) {
        P[GenotypeSource::VCF_CALL] = calculateGenotypeFrequencies(gt);
        this->genotypes[GenotypeSource::VCF_CALL] = gt;
    }

    inline void setFilter(Filter f) { this->filter = f; }

    inline std::string getChromosome() { return this->chrom; }
    inline std::string getRef() { return this->ref; }
    inline std::string getAlt() { return this->alt; }

    inline int getPosition() { return this->pos; }

    inline std::vector<GenotypeSource> getAllGenotypes(){
        std::vector<GenotypeSource> all;
        for(auto const& gt : this->genotypes)
            all.push_back(gt.first);
        return all;
    }
    //inline VectorXd getGenotype(Genotype gt) { return genotypes[gt]; }
    inline VectorXd* getGenotype(GenotypeSource gt) { return &genotypes[gt]; }

    inline Vector3d* getP(GenotypeSource gt) { return &P[gt]; }

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

    inline void shrink(){
        genotypes.clear();
        P.clear();
        isShrunk = true;
    }

};

inline bool variantCompare(Variant lhs, Variant rhs) { return lhs < rhs; }

struct VariantSet{
private:
    std::vector<Variant> variants;
    std::vector<double> pval;
    int nvalid = 0;
    Interval *interval;
    bool hasInterval = false;
public:
    VariantSet(Variant& v) { variants.push_back(v); if(v.isValid()) nvalid++; }
    VariantSet() { }

    inline void setInterval(Interval * inv) { interval = inv; hasInterval = true;}
    inline bool isIn(Variant &variant) { return interval->isIn(variant.getChromosome(), variant.getPosition()); }

    inline void addVariant(Variant &variant) {
        variants.push_back(variant);
        if(variant.isValid()) nvalid++;
    }

    inline std::vector<Variant>* getVariants() { return &variants; }
    inline std::string getChromosome() { return (variants.size() < 1) ? "na" : variants[0].getChromosome(); }
    inline int getMinPos() {
        for(size_t i = 0; i < variants.size(); i++)
            if(variants[i].isValid())
                return variants[i].getPosition();

        return -1;
    }
    inline int getMaxPos() {
        for(size_t i = variants.size()-1; i <= 0; i++)
            if(variants[i].isValid())
                return variants[i].getPosition();

        return -1;
    }
    inline double getMidPos() { return (variants.size() < 1) ? -1 :
                                                            (variants[0].getPosition() + variants.back().getPosition())/2.0; }

    inline void setPval(int i, double p) { pval[i] = p; }

    inline void addPval(double p) { pval.push_back(p); }
    inline double getPval(size_t i) { return pval[i]; }

    inline int nPvals() { return static_cast<int>(this->pval.size()); }
    inline int size(){ return static_cast<int>(variants.size()); }
    inline int validSize(){ return nvalid; }
    inline bool isValid(){ return nvalid > 0; }

    inline MatrixXd getX(GenotypeSource gt){
        std::vector<VectorXd*> x;
        for (size_t i = 0; i < variants.size(); i++)
            if(variants[i].isValid())
                x.push_back(variants[i].getGenotype(gt));

        if(x.size() < 1)
            return MatrixXd();

        MatrixXd X(x[0]->rows(), x.size());
        for (size_t i = 0; i < x.size(); i++)
            X.col(static_cast<int>(i)) = *x[i];

        return X;
    }

    inline MatrixXd getP(GenotypeSource gt){
        std::vector<Vector3d*> p;
        for (size_t i = 0; i < variants.size(); i++)
            if(variants[i].isValid())
                p.push_back(variants[i].getP(gt));

        if(p.size() < 1)
            return MatrixXd();

        MatrixXd P(p.size(), 3);
        for (size_t i = 0; i < p.size(); i++)
            P.row(static_cast<int>(i)) = *p[i];

        return P;
    }

    inline void shrink(){
        for (size_t i = 0; i < variants.size(); i++)
            variants[i].shrink();
    }

    inline std::string toString(std::string test, int pvalIndex, int setID){

        std::string str = "";

        if(variants.size() < 1 || pvalIndex >= nPvals())
            return  str;

        for(size_t i = 0; i < variants.size(); i++){
            str += variants[i].toString();
            str += "\t" + std::to_string(pval[pvalIndex]);
            str += "\t" + test;

            if(variants.size() > 1){
                if(hasInterval)
                    str += "\t" + interval->id;
                else
                    str += "\t" + std::to_string(setID);
            }

            if(i != variants.size() - 1)
                str+= "\n";
        }

        return str;
    }
};
