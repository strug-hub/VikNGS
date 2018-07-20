#include "tabledisplaywindow.h"
#include "ui_tabledisplaywindow.h"


void TableDisplayWindow::addSampleInfo(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, Variant& variant){

    titles.append("index");
    titles.append("Cohort");
    titles.append("Group ID");

    if(simulation)
        titles.append("Read Depth");

    QVector<QTableWidgetItem*> index(nrow);
    QVector<QTableWidgetItem*> caseControl(nrow);
    QVector<QTableWidgetItem*> groupID(nrow);
    QVector<QTableWidgetItem*> readDepth(nrow);

    for (int i = 0; i < nrow ; i++){
        index[i] =  new QTableWidgetItem(QString::number(i));
        caseControl[i] =  new QTableWidgetItem((y[i] == 1? "case":"control"));
        groupID[i] =  new QTableWidgetItem(QString::number(g[i]));

        if(simulation)
            readDepth[i] = new QTableWidgetItem(QString::number(variant.readDepths[i]));
    }

    table.push_back(index);
    table.push_back(caseControl);
    table.push_back(groupID);
    if(simulation)
        table.push_back(readDepth);

    return;
}

void TableDisplayWindow::addSampleGenotypes(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, Variant& variant){

    if(simulation)
        titles.append("True GT");

    titles.append("Called GT");
    titles.append("Expected GT");

    QVector<QTableWidgetItem*> trueGT(nrow);
    QVector<QTableWidgetItem*> calledGT(nrow);
    QVector<QTableWidgetItem*> expectedGT(nrow);

    QColor gtColour = QColor(194, 214, 211);

    for (int i = 0; i < nrow ; i++){

        if(simulation){
            QTableWidgetItem* tcell = new QTableWidgetItem(QString::number(variant.trueGenotype[i]));
            gtColour.setAlpha( std::min(255.0, (variant.trueGenotype[i]/2)*255 ) );
            tcell->setBackgroundColor(gtColour);
            trueGT[i]=tcell;
        }

        QTableWidgetItem* ccell = new QTableWidgetItem(QString::number(variant.genotypeCalls[i]));
        gtColour.setAlpha( std::min(255.0, (variant.genotypeCalls[i]/2)*255 ) );
        ccell->setBackgroundColor(gtColour);
        calledGT[i]=ccell;

        QTableWidgetItem* ecell = new QTableWidgetItem(QString::number(variant.expectedGenotype[i]));
        gtColour.setAlpha( std::min(255.0, (variant.expectedGenotype[i]/2)*255 ) );
        ecell->setBackgroundColor(gtColour);
        expectedGT[i]=ecell;
    }

    if(simulation)
        table.push_back(trueGT);
    table.push_back(calledGT);
    table.push_back(expectedGT);

    return;
}

void TableDisplayWindow::addVCFData(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, Variant& variant){

    if(simulation)
       return;

    titles.append("INFO");
    titles.append("FORMAT");
    titles.append("COLUMN USED");

    QVector<QTableWidgetItem*> format(nrow);
    QVector<QTableWidgetItem*> info(nrow);
    QVector<QTableWidgetItem*> used(nrow);


    for (int i = 0; i < nrow ; i++){
        format[i] = new QTableWidgetItem(QString::fromStdString(variant.format));
        info[i] = new QTableWidgetItem(QString::fromStdString(variant.vcfCalls[i]));
        used[i] = new QTableWidgetItem(QString::fromStdString(variant.columnUsed[i]));

    }

    table.push_back(info);
    table.push_back(format);
    table.push_back(used);

    return;
}

void TableDisplayWindow::addInfo(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table){

    titles.append("index");
    titles.append("CHROM");
    titles.append("POS");
    titles.append("REF");
    titles.append("ALT");

    //    CheckTree* info;
    //    info->initalize("Variant Info");
    //    info->addChild("CHROM", 0);
     //   info->addChild("POS", 1);
      //  info->addChild("REF", 2);
     //   info->addChild("ALT", 3);
     //   variantCheckTree.push_back(info);


    QVector<QTableWidgetItem*> index(nrow);
    QVector<QTableWidgetItem*> chr(nrow);
    QVector<QTableWidgetItem*> pos(nrow);
    QVector<QTableWidgetItem*> ref(nrow);
    QVector<QTableWidgetItem*> alt(nrow);

    for (int i = 0; i < nrow ; i++){
        index[i] = new QTableWidgetItem(QString::number(i));
        chr[i] = new QTableWidgetItem(QString::fromStdString(variants->at(i).chr));
        pos[i] = new QTableWidgetItem(QString::number(variants->at(i).pos));
        ref[i] = new QTableWidgetItem(QString::fromStdString(variants->at(i).ref));
        alt[i] = new QTableWidgetItem(QString::fromStdString(variants->at(i).alt));
    }

    table.push_back(index);
    table.push_back(chr);
    table.push_back(pos);
    table.push_back(ref);
    table.push_back(alt);

    return;
}

void TableDisplayWindow::addAlleleProbability(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, Variant& variant){

    titles.append("p_0");
    titles.append("p_1");
    titles.append("p_2");

    QVector<QTableWidgetItem*> p0(nrow);
    QVector<QTableWidgetItem*> p1(nrow);
    QVector<QTableWidgetItem*> p2(nrow);

    for (int i = 0; i < nrow ; i++){
        p0[i] = new QTableWidgetItem(QString::number(variant.P[0]));
        p1[i] = new QTableWidgetItem(QString::number(variant.P[1]));
        p2[i] = new QTableWidgetItem(QString::number(variant.P[2]));
    }

    table.push_back(p0);
    table.push_back(p1);
    table.push_back(p2);

    return;
}

void TableDisplayWindow::addPvals(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table){

    QColor pvalColour = QColor(234, 183, 53, 255);

    for (int j = 0; j < variants->at(0).nPvals(); j++){
        titles.append(QString::fromStdString(variants->at(0).getPvalSourceShort(j)) + " pval");

        QVector<QTableWidgetItem*> pvals(nrow);
        for (int i = 0; i < nrow ; i++){
            double pval = variants->at(i).getPval(j);
            pvalColour.setAlpha( std::min(255.0, (-log(pval+1e-14)/10)*255 ));
            QTableWidgetItem* pcell = new QTableWidgetItem(QString::number(pval));
            pcell->setBackgroundColor(pvalColour);
            pvals[i] = pcell;
        }

        table.push_back(pvals);
    }

    return;
}

double TableDisplayWindow::calculateMaf(VectorXd& gt, bool caseOnly, bool controlOnly){

    double maf;

    double denom = 0;

    for (int i = 0; i < gt.rows() ; i++){

        if(std::isnan(gt[i]))
            continue;
        else if(caseOnly && y[i] == 0)
            continue;
        else if(controlOnly && y[i] == 1)
            continue;

        maf += gt[i];
        denom+=2;
    }

    return maf/std::max(1.0,denom);
}


void TableDisplayWindow::addMafs(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, QString gt, QString caseControl){

    if(gt=="true" && !simulation)
        return;

    QString title = "MAF (";
    if(gt == "true")
        title.append("true");
    else if(gt == "calls")
        title.append("call");
    else if(gt == "expected")
        title.append("exp");

    if(caseControl == "case")
        title.append(", case)");
    else if(gt == "control")
        title.append(", control)");
    else
        title.append(")");

    titles.append(title);

    QVector<QTableWidgetItem*> mafs(nrow);
    QColor mafColour = QColor(194, 214, 211);

    for (int i = 0; i < nrow ; i++){

        double maf;

        if(gt == "true")
            maf = calculateMaf(variants->at(i).trueGenotype, caseControl=="case", caseControl=="control");
        else if (gt == "calls")
            maf = calculateMaf(variants->at(i).trueGenotype, caseControl=="case", caseControl=="control");
        else if(gt == "expected")
            maf = calculateMaf(variants->at(i).expectedGenotype, caseControl=="case", caseControl=="control");

        mafColour.setAlpha(std::min(255.0, maf*255*2));
        QTableWidgetItem* mcell = new QTableWidgetItem(QString::number(maf));
        mcell->setBackgroundColor(mafColour);
        mafs[i] = mcell;

    }

    table.push_back(mafs);

    return;
}


void TableDisplayWindow::addGtFreq(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, QString gt, QString caseControl){

    if(gt=="true" && !simulation)
        return;

    QString suffix = "Freq (";
    if(gt == "true")
        suffix.append("true");
    else if(gt == "calls")
        suffix.append("call");

    if(caseControl == "case")
        suffix.append(", case)");
    else if(gt == "control")
        suffix.append(", control)");
    else
        suffix.append(")");

    titles.append(QString("0_").append(suffix));
    titles.append(QString("1_").append(suffix));
    titles.append(QString("2_").append(suffix));

    QVector<QTableWidgetItem*> freq0(nrow);
    QVector<QTableWidgetItem*> freq1(nrow);
    QVector<QTableWidgetItem*> freq2(nrow);

    for (int i = 0; i < nrow ; i++){

        if(caseControl == "case" && y[i]==0)
            continue;
        if(caseControl == "control" && y[i]==1)
            continue;

        int _0 = 0;
        int _1 = 0;
        int _2 = 0;
        double genotype = -1;

        if(gt == "true"){
            for (int j = 0; j < variants->at(i).trueGenotype.rows() ; j++){
                genotype = variants->at(i).trueGenotype[j];
                if(genotype == 0) _0++;
                else if (genotype == 1) _1++;
                else _2++;
               }
        }

        else if(gt == "calls"){
            for (int j = 0; j < variants->at(i).genotypeCalls.rows() ; j++){
                genotype = variants->at(i).genotypeCalls[j];
                if(genotype == 0) _0++;
                else if (genotype == 1) _1++;
                else _2++;
            }
        }

        freq0[i] = new QTableWidgetItem(QString::number(_0));
        freq1[i] = new QTableWidgetItem(QString::number(_1));
        freq2[i] = new QTableWidgetItem(QString::number(_2));
    }

    table.push_back(freq0);
    table.push_back(freq1);
    table.push_back(freq2);

    return;
}


