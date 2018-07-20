#include "tabledisplaywindow.h"
#include "ui_tabledisplaywindow.h"

TableDisplayWindow::TableDisplayWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::TableDisplayWindow)
{
    ui->setupUi(this);

    setWindowIcon(QIcon(":icon.svg"));

}

TableDisplayWindow::~TableDisplayWindow()
{
    delete ui;
}

void TableDisplayWindow::initialize(QString title, std::vector<Variant>& variants, SimulationRequest& simRequest, int index){

    simulation = true;

    this->setWindowTitle("Table - " + title);
    this->variants = &variants;
    this->simRequest = &simRequest;
    this->g = simRequest.getGroups(index);
    this->y = simRequest.getCaseControlStatus(index);
    buildVariantTable();
}


void TableDisplayWindow::initialize(QString title, std::vector<Variant>& variants, TestInput& input, int index){

    simulation = false;
    this->setWindowTitle("Table - " + title);
    this->variants = &variants;
    this->input = &input;
    //todo
    //this->g = request.getGroups(index);
    //this->y = request.getCaseControlStatus(index);
    buildVariantTable();
}


void TableDisplayWindow::fillGenotypeTable(int variantIndex){

    ui->table_genoTbl->setColumnCount(0);

    Variant v = variants->at(variantIndex);
    int nrow = v.trueGenotype.rows();

    QStringList titles;

    titles.append("index");
    titles.append("Cohort");
    titles.append("Group ID");
    titles.append("Read Depth");

    QVector<QString> sampleID(nrow);
    QVector<QString> caseControl(nrow);
    QVector<QString> groupID(nrow);
    QVector<int> readDepth(nrow);

    titles.append("True GT");
    titles.append("Called GT");
    titles.append("Expected GT");

    QVector<QVector<double>> genotypes(nrow);

    QColor gtColour = QColor(194, 214, 211);

    for (int i = 0; i < nrow ; i++){

        sampleID[i] = QString::number(i);
        caseControl[i] = (y[i] == 1? "case":"control");
        groupID[i] = QString::number(g[i]);

        QVector<double> gt;

        gt.push_back(v.trueGenotype[i]);
        gt.push_back(v.genotypeCalls[i]);
        gt.push_back(v.expectedGenotype[i]);
        genotypes[i] = gt;
        readDepth[i] = v.readDepths[i];

    }

    int ncol = titles.size();
    ui->table_genoTbl->setColumnCount(ncol);
    ui->table_genoTbl->setRowCount(nrow);
    for (int i = 0; i < nrow ; i++){

        ui->table_genoTbl->setItem(i, 0, new QTableWidgetItem(sampleID[i]));
        ui->table_genoTbl->setItem(i, 1, new QTableWidgetItem(caseControl[i]));
        ui->table_genoTbl->setItem(i, 2, new QTableWidgetItem(groupID[i]));

        ui->table_genoTbl->setItem(i, 3, new QTableWidgetItem(QString::number(readDepth[i])));

        int ngenotypes = genotypes[i].size();
        for (int j = 0; j < ngenotypes; j++){
            gtColour.setAlpha( std::min(255.0, (genotypes[i][j]/2)*255 ) );
            QTableWidgetItem* gtcell = new QTableWidgetItem(QString::number(genotypes[i][j]));
            gtcell->setBackgroundColor(gtColour);
            ui->table_genoTbl->setItem(i, 4+j, gtcell);
        }
    }

    ui->table_genoTbl->setHorizontalHeaderLabels(titles);
}


void TableDisplayWindow::drawVariantCheckTree(){

    for(CheckTree* tree : variantCheckTree){
        QTreeWidgetItem* branch = tree->draw();
        ui->table_checkTre->addTopLevelItem(branch);
        connect(ui->table_checkTre, SIGNAL(itemSelectionChanged()),
                tree, SLOT(updateState()));
    }
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

void TableDisplayWindow::addMafs(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, QString gt, QString caseControl){

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

    int ncase = y.sum();
    int ncontrol = (1-y.array()).sum();

    for (int i = 0; i < nrow ; i++){

        double maf;
        if(caseControl == "case" && gt == "true")
            maf = (variants->at(i).trueGenotype.array() * y.array()).sum() / (2*ncase);
        else if(caseControl == "case" && gt == "expected")
            maf = (variants->at(i).expectedGenotype.array() * y.array()).sum() / (2*ncase);
        else if(caseControl == "case" && gt == "calls")
            maf = (variants->at(i).genotypeCalls.array() * y.array()).sum() / (2*ncase);
        else if(caseControl == "control" && gt == "true")
            maf = (variants->at(i).trueGenotype.array() * (1-y.array())).sum() / (2*ncontrol);
        else if(caseControl == "control" && gt == "calls")
            maf = (variants->at(i).genotypeCalls.array() * (1-y.array())).sum() / (2*ncontrol);
        else if(caseControl == "control" && gt == "expected")
            maf = (variants->at(i).expectedGenotype.array() * (1-y.array())).sum() / (2*ncontrol);
        else if(gt == "true")
            maf = variants->at(i).trueGenotype.mean()/2;
        else if (gt == "expected")
            maf = variants->at(i).expectedGenotype.mean()/2;
        else if (gt == "calls")
            maf = variants->at(i).genotypeCalls.mean()/2;

        mafColour.setAlpha(std::min(255.0, maf*255*2));
        QTableWidgetItem* mcell = new QTableWidgetItem(QString::number(maf));
        mcell->setBackgroundColor(mafColour);
        mafs[i] = mcell;

    }

    table.push_back(mafs);

    return;
}


void TableDisplayWindow::addGtFreq(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, QString gt, QString caseControl){

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


void TableDisplayWindow::buildVariantTable(){

    ui->table_variantTbl->clearContents();
    int nrow = variants->size();

    QStringList titles;
    QVector<QVector<QTableWidgetItem*>> table;

    addInfo(nrow, titles, table);
    addPvals(nrow, titles, table);
    addMafs(nrow, titles, table, "true");
    addMafs(nrow, titles, table, "calls");
    addMafs(nrow, titles, table, "expected");
    addMafs(nrow, titles, table, "true", "case");
    addMafs(nrow, titles, table, "true", "control");
    addMafs(nrow, titles, table, "calls", "case");
    addMafs(nrow, titles, table, "calls", "control");
    addMafs(nrow, titles, table, "expected", "case");
    addMafs(nrow, titles, table, "expected", "control");

    addGtFreq(nrow, titles, table, "true");
    addGtFreq(nrow, titles, table, "calls");

    int ncol = titles.size();
    ui->table_variantTbl->setColumnCount(ncol);
    ui->table_variantTbl->setRowCount(nrow);
    ui->table_variantTbl->setHorizontalHeaderLabels(titles);

    for(int i = 0; i < table.size(); i++)
        for(int j = 0; j < table[i].size(); j++)
            ui->table_variantTbl->setItem(j,i,table[i][j]);


    //  drawVariantCheckTree();*/
}

void TableDisplayWindow::displayCalls(int index){
    QString ref = QString::fromStdString(variants->at(selectedVariantIndex).ref);
    QString alt = QString::fromStdString(variants->at(selectedVariantIndex).alt);
    QString error1;
    QString error2;

    if((ref == "A" || ref=="T") && (alt == "A" || alt == "T")){
        error1 = "G";
        error2 = "C";
    }
    else if((ref == "A" || ref=="C") && (alt == "A" || alt == "C")){
        error1 = "G";
        error2 = "T";
    }
    else if((ref == "A" || ref=="G") && (alt == "A" || alt == "G")){
        error1 = "T";
        error2 = "C";
    }
    else if((ref == "T" || ref=="C") && (alt == "T" || alt == "C")){
        error1 = "A";
        error2 = "G";
    }
    else if((ref == "T" || ref=="G") && (alt == "T" || alt == "G")){
        error1 = "A";
        error2 = "C";
    }
    else if((ref == "C" || ref=="G") && (alt == "C" || alt == "G")){
        error1 = "A";
        error2 = "T";
    }

    //[ref, alt, error1, error2]
    std::vector<int> baseCalls = variants->at(selectedVariantIndex).baseCalls[index];

    QString calls = "";
    for(int i = 0; i < baseCalls[0]; i++)
        calls.append(basePair(ref));
    for(int i = 0; i < baseCalls[1]; i++)
        calls.append(basePair(alt));
    for(int i = 0; i < baseCalls[2]; i++)
        calls.append(basePair(error1));
    for(int i = 0; i < baseCalls[3]; i++)
        calls.append(basePair(error2));

    ui->table_callsTxt->setText(calls);

}

void TableDisplayWindow::on_table_variantTbl_cellClicked(int row, int column){

    int index = ui->table_variantTbl->item(row,0)->text().toInt();
    fillGenotypeTable(index);

    QString ref = QString::fromStdString(variants->at(index).ref);
    QString alt = QString::fromStdString(variants->at(index).alt);
    ui->table_refLbl->setText(basePair(ref));
    ui->table_altLbl->setText(basePair(alt));
    selectedVariantIndex = index;
}

void TableDisplayWindow::on_table_genoTbl_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn)
{
    if(currentRow >= 0){
        int index = ui->table_genoTbl->item(currentRow,0)->text().toInt();
        displayCalls(index);
    }
}
