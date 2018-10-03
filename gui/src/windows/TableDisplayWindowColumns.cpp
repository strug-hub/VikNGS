#include "TableDisplayWindow.h"
#include "ui_tabledisplaywindow.h"


void TableDisplayWindow::addSampleInfo(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table){
    titles.append("index");
    titles.append("Cohort");
    titles.append("Group ID");

    QVector<QTableWidgetItem*> index(nrow);
    QVector<QTableWidgetItem*> caseControl(nrow);
    QVector<QTableWidgetItem*> groupID(nrow);
    QVector<QTableWidgetItem*> readDepth(nrow);

    for (int i = 0; i < nrow ; i++){
        index[i] =  new QTableWidgetItem(QString::number(i));
        caseControl[i] =  new QTableWidgetItem((abs(y[i]) < 1e-4 ? "control":"case"));
        groupID[i] =  new QTableWidgetItem(QString::number(g[i]));
    }

    table.push_back(index);
    table.push_back(caseControl);
    table.push_back(groupID);
}

void TableDisplayWindow::addSampleGenotypes(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, Variant* variant){

    QColor gtColour = QColor(194, 214, 211);

    VectorXd* X;
    for(Genotype gt : variant->getAllGenotypes()){
        titles.append(QString::fromStdString(genotypeToString(gt) + " GT"));
        X = variant->getGenotype(gt);
        QVector<QTableWidgetItem*> gtValues(nrow);

        for(int i = 0; i < nrow; i++){

            QTableWidgetItem* gtcell = new QTableWidgetItem(QString::number(X->coeff(i)));
            gtColour.setAlpha( std::min(255.0, (X->coeff(i)/2)*255 ) );
            gtcell->setBackgroundColor(gtColour);
            gtValues[i]=gtcell;
        }
        table.push_back(gtValues);
    }
}

void TableDisplayWindow::addInfo(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table){

    titles.append("index");
    titles.append("CHROM");
    titles.append("POS");
    titles.append("REF");
    titles.append("ALT");

    QVector<QTableWidgetItem*> index(nrow);
    QVector<QTableWidgetItem*> chr(nrow);
    QVector<QTableWidgetItem*> pos(nrow);
    QVector<QTableWidgetItem*> ref(nrow);
    QVector<QTableWidgetItem*> alt(nrow);

    int i = 0;
    for(size_t n = 0; n < variants->size(); n++){
        for(size_t m = 0; m < variants->at(n).size(); m++){
            index[i] = new QTableWidgetItem(QString::number(i));
            chr[i] = new QTableWidgetItem(QString::fromStdString(variants->at(n).getVariants()->at(m).getChromosome()));
            pos[i] = new QTableWidgetItem(QString::number(variants->at(n).getVariants()->at(m).getPosition()));
            ref[i] = new QTableWidgetItem(QString::fromStdString(variants->at(n).getVariants()->at(m).getRef()));
            alt[i] = new QTableWidgetItem(QString::fromStdString(variants->at(n).getVariants()->at(m).getAlt()));
            i++;
        }
    }

    table.push_back(index);
    table.push_back(chr);
    table.push_back(pos);
    table.push_back(ref);
    table.push_back(alt);
}

void TableDisplayWindow::addPvals(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table){

    QColor pvalColour = QColor(234, 183, 53, 255);

    for (size_t j = 0; j < tests->size(); j++){

        bool contains = false;
        for(int idx : testIndexes)
            if (idx == static_cast<int>(j))
                contains = true;

        if(!contains)
            continue;

        titles.append(QString::fromStdString(tests->at(j).toShortString()) + " pval");
        QVector<QTableWidgetItem*> pvals(nrow);
        int i = 0;

        for(size_t n = 0; n < variants->size(); n++){
            for(size_t m = 0; m < variants->at(n).size(); m++){
                double pval = variants->at(n).getPval(j);
                pvalColour.setAlpha( std::min(255.0, (-log(pval+1e-14)/10)*255.0 ));
                QTableWidgetItem* pcell = new QTableWidgetItem(QString::number(pval));
                pcell->setBackgroundColor(pvalColour);
                pvals[i] = pcell;
                i++;
            }
        }

        table.push_back(pvals);
    }
}

double TableDisplayWindow::calculateMaf(VectorXd* gt, bool caseOnly, bool controlOnly){

    double maf = 0;
    double denom = 0;

    int sampleSize = nsamples;
    if (sampleSize < 0) sampleSize = gt->rows();
    if(sampleSize > gt->rows()) sampleSize = gt->rows();

    for (int i = 0; i < sampleSize; i++){

        if(std::isnan(gt->coeff(i)))
            continue;
        else if(caseOnly && abs(y[i]) < 1e-4)
            continue;
        else if(controlOnly && abs(y[i]-1) < 1e-4)
            continue;

        maf += gt->coeff(i);
        denom+=2;
    }

    return maf/std::max(1.0,denom);
}

void TableDisplayWindow::addMafsCaseControl(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, bool useCases){

    QColor mafColour = QColor(194, 214, 211);
    std::string caseControl = useCases ? " case MAF" : " control MAF";
    for(size_t j = 0; j < testIndexes.size(); j++){
        Genotype genotype = tests->at(testIndexes[j]).getGenotype();
        titles.append(QString::fromStdString(genotypeToString(genotype) + caseControl));

        QVector<QTableWidgetItem*> mafs(nrow);

        int i = 0;
        for(size_t n = 0; n < variants->size(); n++){
            for(size_t m = 0; m < variants->at(n).size(); m++){

                double maf = calculateMaf(variants->at(n).getVariants()->at(m).getGenotype(genotype), useCases, !useCases);
                mafColour.setAlpha(std::min(255.0, maf*255.0*2.0));
                QTableWidgetItem* mcell = new QTableWidgetItem(QString::number(maf));
                mcell->setBackgroundColor(mafColour);
                mafs[i] = mcell;
                i++;
            }
        }
        table.push_back(mafs);
    }
}

void TableDisplayWindow::addMafs(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table){

    QColor mafColour = QColor(194, 214, 211);
    for(size_t j = 0; j < tests->size(); j++){
        Genotype genotype = tests->at(j).getGenotype();
        titles.append(QString::fromStdString(genotypeToString(genotype) + " MAF"));

        QVector<QTableWidgetItem*> mafs(nrow);

        int i = 0;
        for(size_t n = 0; n < variants->size(); n++){
            for(size_t m = 0; m < variants->at(n).size(); m++){

                double maf = calculateMaf(variants->at(n).getVariants()->at(m).getGenotype(genotype), false, false);
                mafColour.setAlpha(std::min(255.0, maf*255.0*2.0));
                QTableWidgetItem* mcell = new QTableWidgetItem(QString::number(maf));
                mcell->setBackgroundColor(mafColour);
                mafs[i] = mcell;
                i++;
            }
        }
        table.push_back(mafs);
    }

}
