#pragma once
#include <QColor>
#include <QVector>
#include <QString>
#include "../src/Variant.h"

class Chromosome{
private:
    QVector<Variant*> variants;
    QVector<VariantSet*> collapse;
    QColor colour;
    QString name;
    double maxPos;
    double minPos;
    double offset;

public:
    Chromosome(){
        name="unknown";
        maxPos = -1;
        minPos = 2100000000;
    }

    Chromosome(QString name, Variant* first,  VariantSet* set, QColor c){
        this->name = name;
        variants.push_back(first);
        collapse.push_back(set);
        maxPos = first->getPosition();
        minPos = first->getPosition();
        colour = c;
    }

    void addVariant(Variant* toAdd, VariantSet* set){
        variants.push_back(toAdd);
        collapse.push_back(set);
        if(maxPos < variants.back()->getPosition())
            maxPos = variants.back()->getPosition();
        if(minPos > variants.back()->getPosition())
            minPos = variants.back()->getPosition();
    }

    QVector<double> getRelativePositions(double offset){
        QVector<double> pos;
        double min = getMinPos();
        for(int i = 0; i < variants.size(); i++)
            pos.push_back(offset + variants[i]->getPosition() - min);

        return pos;
    }

    QVector<double> getPositions(){
        QVector<double> pos;
        for(int i = 0; i < variants.size(); i++)
            pos.push_back(variants[i]->getPosition());

        return pos;
    }

    inline double getPosition(int index){ return variants[index]->getPosition(); }

    QVector<double> getPvals(int index=0){
        QVector<double> pvals; pvals.reserve(collapse.size());
        for(VariantSet* vs : collapse)
            pvals.push_back(-log10(vs->getPval(index)));
        return pvals;
    }

    inline double getPval(int index){
        return -log10(collapse[index]->getPval(0));
    }

    Variant* getVariant(int i){return variants[i];}
    void setColour(QColor c){ this->colour = c; }
    double getOffset(){ return this->offset; }
    void setOffset(double offset){ this->offset = offset; }
    QColor getColour(){ return this->colour; }
    int size(){ return variants.size(); }
    double getSpan(){ return maxPos - minPos; }
    double getMaxPos(){return maxPos; }
    double getMinPos(){ return minPos; }
    QString getName(){ return name; }
    void setName(QString n){ this->name = n; }

};

inline int extractNumber(const QString &c)
{
    QString number = "";

    for(int i=0; i < c.size(); i++)
        if(c[i].isDigit())
            number += c[i];

    if(number == "")
        return -1;
    return number.toInt();

}

inline bool chrCompare(const QString &c1, const QString &c2)
{
    int i1 = extractNumber(c1);
    int i2 = extractNumber(c2);

    if(i1 > 0 && i2 > 0)
        return i1 < i2;

    return c1.toLower() < c2.toLower();
}

