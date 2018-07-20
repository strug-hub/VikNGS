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

void TableDisplayWindow::initialize(QString title, std::vector<Variant>& variants, SimulationRequest& request, int index){

    this->setWindowTitle("Table - " + title);
    this->variants = &variants;
    this->request = &request;
    this->g = request.getGroups(index);
    this->y = request.getCaseControlStatus(index);
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


void TableDisplayWindow::buildVariantTable(){

    ui->table_variantTbl->clearContents();
    int nrow = variants->size();

    QColor pvalColour = QColor(234, 183, 53, 255);
    QColor mafColour = QColor(194, 214, 211);

    QStringList titles;

//    CheckTree* info;
//    info->initalize("Variant Info");
//    info->addChild("CHROM", 0);
 //   info->addChild("POS", 1);
  //  info->addChild("REF", 2);
 //   info->addChild("ALT", 3);
 //   variantCheckTree.push_back(info);

    titles.append("index");
    QVector<QString> index(nrow);


    titles.append("CHROM");
    titles.append("POS");
    titles.append("REF");
    titles.append("ALT");
    QVector<QString> chr(nrow);
    QVector<QString> pos(nrow);
    QVector<QString> ref(nrow);
    QVector<QString> alt(nrow);

    QVector<QVector<double>> pvals(nrow);
    QVector<QVector<double>> mafs(nrow);


    QVector<QVector<QString>> freq(nrow);

    QVector<double> miscall(nrow);

    titles.append("Miscall");

    for (int i = 0; i < variants->at(0).nPvals(); i++)
        titles.append(QString::fromStdString(variants->at(0).getPvalSourceShort(i)) + " pval");



    titles.append("MAF (true)");
    titles.append("MAF (sample)");
    titles.append("MAF (GT calls)");
    titles.append("MAF (expected GT)");

    titles.append("MAF (true, case)");
    titles.append("MAF (true, control)");
    titles.append("MAF (call, case)");
    titles.append("MAF (call, control)");
    titles.append("MAF (exp, case)");
    titles.append("MAF (exp, control)");

    titles.append("GT Freq 0");
    titles.append("GT Freq 1");
    titles.append("GT Freq 2");


    for (int i = 0; i < nrow ; i++){
        index[i] = QString::number(i);

        chr[i] = (QString::fromStdString(variants->at(i).chr));
        pos[i] = (QString::number(variants->at(i).pos));
        ref[i] = (QString::fromStdString(variants->at(i).ref));
        alt[i] =(QString::fromStdString(variants->at(i).alt));
        QVector<double> pval;
        for (int j = 0; j < variants->at(i).nPvals(); j++)
            pval.push_back(variants->at(i).getPval(j));
        pvals[i] = pval;

        QVector<double> maf;

        maf.push_back(variants->at(i).trueMaf);
        maf.push_back(variants->at(i).trueGenotype.mean()/2);
        maf.push_back(variants->at(i).genotypeCalls.mean()/2);
        maf.push_back(variants->at(i).expectedGenotype.mean()/2);

        int ncase = y.sum();
        int ncontrol = (1-y.array()).sum();

        maf.push_back(((variants->at(i).trueGenotype.array() * y.array()).sum()/2) / ncase);
        maf.push_back(((variants->at(i).trueGenotype.array() * (1-y.array())).sum()/2) / ncontrol);
        maf.push_back(((variants->at(i).genotypeCalls.array() * y.array()).sum()/2) / ncase);
        maf.push_back(((variants->at(i).genotypeCalls.array() * (1-y.array())).sum()/2) / ncontrol);
        maf.push_back(((variants->at(i).expectedGenotype.array() * y.array()).sum()/2 / ncase));
        maf.push_back(((variants->at(i).expectedGenotype.array() * (1-y.array())).sum()/2 / ncontrol));


        mafs[i] = maf;

        miscall[i] =(variants->at(i).genotypeCalls.array() - variants->at(i).trueGenotype.array()).abs().sum();

        QVector<QString> fq;
        int freq0 = 0;
        int freq1 = 0;
        int freq2 = 0;

        int freq0t = 0;
        int freq1t = 0;
        int freq2t = 0;

        for(int j = 0; j < variants->at(i).genotypeCalls.rows(); j++){
            freq0t +=(variants->at(i).trueGenotype[j]==0? 1:0);
            freq1t +=(variants->at(i).trueGenotype[j]==1? 1:0);
            freq2t +=(variants->at(i).trueGenotype[j]==2? 1:0);

            freq0+=(variants->at(i).genotypeCalls[j]==0? 1:0);
            freq1+=(variants->at(i).genotypeCalls[j]==1? 1:0);
            freq2+=(variants->at(i).genotypeCalls[j]==2? 1:0);


        }

        fq.push_back(QString::number(freq0t).append("/").append(QString::number(freq0)));
        fq.push_back(QString::number(freq1t).append("/").append(QString::number(freq1)));
        fq.push_back(QString::number(freq2t).append("/").append(QString::number(freq2)));

        freq[i] = fq;
    }

    int ncol = titles.size();
    ui->table_variantTbl->setColumnCount(ncol);
    ui->table_variantTbl->setRowCount(nrow);
    for (int i = 0; i < nrow ; i++){
        ui->table_variantTbl->setItem(i, 0, new QTableWidgetItem(index[i]));

        ui->table_variantTbl->setItem(i, 1, new QTableWidgetItem(chr[i]));
        ui->table_variantTbl->setItem(i, 2, new QTableWidgetItem(pos[i]));
        ui->table_variantTbl->setItem(i, 3, new QTableWidgetItem(ref[i]));
        ui->table_variantTbl->setItem(i, 4, new QTableWidgetItem(alt[i]));

        ui->table_variantTbl->setItem(i, 5, new QTableWidgetItem(QString::number(miscall[i])));

        int npval = pvals[i].size();
        for (int j = 0; j < npval; j++){
            pvalColour.setAlpha( std::min(255.0, (-log(pvals[i][j]+1e-14)/10)*255 ));
            QTableWidgetItem* pcell = new QTableWidgetItem(QString::number(pvals[i][j]));
            pcell->setBackgroundColor(pvalColour);
            ui->table_variantTbl->setItem(i, 6+j, pcell);
        }

        int nmafs=mafs[i].size();
        for (int j = 0; j < nmafs; j++){
            mafColour.setAlpha(std::min(255.0, mafs[i][j]*255*2));
            QTableWidgetItem* mcell = new QTableWidgetItem(QString::number(mafs[i][j]));
            mcell->setBackgroundColor(mafColour);
            ui->table_variantTbl->setItem(i, npval+6+j, mcell);
        }

        ui->table_variantTbl->setItem(i, nmafs+npval+6, new QTableWidgetItem(freq[i][0]));
        ui->table_variantTbl->setItem(i, nmafs+npval+6+1, new QTableWidgetItem(freq[i][1]));
        ui->table_variantTbl->setItem(i, nmafs+npval+6+2, new QTableWidgetItem(freq[i][2]));

    }

    ui->table_variantTbl->setHorizontalHeaderLabels(titles);
  //  drawVariantCheckTree();
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
    int index = ui->table_genoTbl->item(currentRow,0)->text().toInt();
    displayCalls(index);
}
