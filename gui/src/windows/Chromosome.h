#pragma once
#include <QColor>
#include <QVector>
#include <QString>
#include "../src/Variant.h"

class Chromosome{
private:
    QVector<VariantSet*> variants;
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

    Chromosome(QString name, VariantSet* first, QColor c){
        this->name = name;
        variants.push_back(first);
        maxPos = first->getMaxPos();
        minPos = first->getMinPos();
        colour = c;
    }

    void addVariant(VariantSet* toAdd){
        variants.push_back(toAdd);
        if(maxPos < toAdd->getMaxPos())
            maxPos = toAdd->getMaxPos();
        if(minPos > toAdd->getMinPos())
            minPos = toAdd->getMinPos();
    }

    QVector<double> getPositions(double offset){

        QVector<double> pos;
        double meanPos;
        for(int i = 0; i < variants.size(); i++){
            meanPos = (variants[i]->getMinPos() + variants[i]->getMaxPos())/2;
            pos.push_back(offset + meanPos);
        }
        return pos;
    }
    QVector<double> getPositions(){ return getPositions(0); }

    QVector<double> getPvals(int index){

        QVector<double> pvals;
        for(int i = 0; i < variants.size(); i++)
            pvals.push_back(-log10(variants[i]->getPval(static_cast<size_t>(index))));

        return pvals;
    }

    VariantSet* getVariant(int i){return variants[i];}
    void setColour(QColor c){ this->colour = c; }
    double getOffset(){ return this->offset; }
    void setOffset(double offset){ this->offset = offset; }
    QColor getColour(){ return this->colour; }
    int size(){ return variants.size(); }
    double getSpan(){ return maxPos - minPos; }
    double getMaxPos(){ return maxPos; }
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

