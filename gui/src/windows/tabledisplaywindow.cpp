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

    this->setWindowTitle("Plotter - " + title);
    this->variants = &variants;
    this->request = &request;
    this->g = request.getGroups(index);
    this->y = request.getCaseControlStatus(index);
    buildVariantTable();
}

void TableDisplayWindow::fillGenotypeTable(int variantIndex){

    ui->table_genoTbl->clearContents();

    Variant v = variants->at(variantIndex);
    int nrow = v.trueGenotype.rows();
    ui->table_genoTbl->setRowCount(nrow);
    ui->table_genoTbl->setColumnCount(4);
    QColor red = Qt::red;

    Vector3d trueFrequency;
    Vector3d calledFrequency;
    double expectedSum;

    for (int i = 0; i < nrow ; i++){

        QString sampleInfo = y[i] == 1? "case":"control";
        sampleInfo.append(" (" + QString::number(g[i]) + ")");
        double trueGT = v.trueGenotype[i];
        trueFrequency[trueGT]++;
        double callGT = v.genotypeCalls[i];
        calledFrequency[callGT]++;
        double expectedGT = v.expectedGenotype[i];
        expectedSum += expectedGT;

        QTableWidgetItem * t = new QTableWidgetItem(QString::number(trueGT));
        QTableWidgetItem * c = new QTableWidgetItem(QString::number(callGT));
        red.setAlpha((std::abs(trueGT-callGT) / 2.0) * 255);
        c->setBackgroundColor(red);

        QTableWidgetItem * e = new QTableWidgetItem(QString::number(expectedGT));
        red.setAlpha((std::abs(trueGT-expectedGT) / 2.0) * 255);
        e->setBackgroundColor(red);


        ui->table_genoTbl->setItem(i, 0,  new QTableWidgetItem(sampleInfo));
        ui->table_genoTbl->setItem(i, 1, t);
        ui->table_genoTbl->setItem(i, 2, c);
        ui->table_genoTbl->setItem(i, 3, e);
    }

    ui->table_genoFreqTbl->clearContents();
    ui->table_genoFreqTbl->setRowCount(3);
    ui->table_genoFreqTbl->setColumnCount(4);

    for (int i = 0; i < 3 ; i++){
        ui->table_genoFreqTbl->setItem(0, i, new QTableWidgetItem(QString::number(trueFrequency[i])));
        ui->table_genoFreqTbl->setItem(1, i, new QTableWidgetItem(QString::number(calledFrequency[i])));
    }

    double averageTrue = (trueFrequency[1] + 2*trueFrequency[2])/(nrow * 1.0);
    double averageCalled = (calledFrequency[1] + 2*calledFrequency[2])/(nrow * 1.0);
    double averageExpected = expectedSum/(nrow * 1.0);
    ui->table_genoFreqTbl->setItem(0, 3, new QTableWidgetItem(QString::number(averageTrue)));
    ui->table_genoFreqTbl->setItem(1, 3, new QTableWidgetItem(QString::number(averageCalled)));
    ui->table_genoFreqTbl->setItem(2, 3, new QTableWidgetItem(QString::number(averageExpected)));
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
    QColor frequencyColour = QColor(194, 214, 211);

    QStringList titles;

    CheckTree* info;
    info->initalize("Variant Info");
    info->addChild("CHROM", 0);
    info->addChild("POS", 1);
    info->addChild("REF", 2);
    info->addChild("ALT", 3);
    variantCheckTree.push_back(info);

    titles.append("CHROM");
    titles.append("POS");
    titles.append("REF");
    titles.append("ALT");
    QVector<QString> chr(nrow);
    QVector<QString> pos(nrow);
    QVector<QString> ref(nrow);
    QVector<QString> alt(nrow);

    QVector<QVector<double>> pvals(nrow);
    QVector<QVector<double>> maf(nrow);

    for (int i = 0; i < variants->at(0).nPvals(); i++)
        titles.append(QString::fromStdString(variants->at(0).getPvalSourceShort(i)) + " pval");

    for (int i = 0; i < nrow ; i++){
        chr[i] = (QString::fromStdString(variants->at(i).chr));
        pos[i] = (QString::number(variants->at(i).pos));
        ref[i] = (QString::fromStdString(variants->at(i).ref));
        alt[i] =(QString::fromStdString(variants->at(i).alt));
        QVector<double> pval;
        for (int j = 0; j < variants->at(i).nPvals(); j++)
            pval.push_back(variants->at(i).getPval(j));
        pvals[i] = pval;
    }

    int ncol = 4 + pvals[0].size();
    ui->table_variantTbl->setColumnCount(ncol);
    ui->table_variantTbl->setRowCount(nrow);
    for (int i = 0; i < nrow ; i++){

        ui->table_variantTbl->setItem(i, 0, new QTableWidgetItem(chr[i]));
        ui->table_variantTbl->setItem(i, 1, new QTableWidgetItem(pos[i]));
        ui->table_variantTbl->setItem(i, 2, new QTableWidgetItem(ref[i]));
        ui->table_variantTbl->setItem(i, 3, new QTableWidgetItem(alt[i]));
        for (int j = 0; j < pvals[i].size(); j++){
            pvalColour.setAlpha( (-log(pvals[i][j]+1e-14)/10)*255 );
            QTableWidgetItem* pcell = new QTableWidgetItem(QString::number(pvals[i][j]));
            pcell->setBackgroundColor(pvalColour);
            ui->table_variantTbl->setItem(i, 4+j, pcell);
        }
    }

    ui->table_variantTbl->setHorizontalHeaderLabels(titles);
    drawVariantCheckTree();
}


void TableDisplayWindow::on_table_variantTbl_cellClicked(int row, int column){
    fillGenotypeTable(row);
}
