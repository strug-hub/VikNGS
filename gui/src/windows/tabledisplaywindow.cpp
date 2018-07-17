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
    fillVariantTable();
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

void TableDisplayWindow::fillVariantTable(){

    ui->table_variantTbl->clearContents();

    int nrow = variants->size();
    ui->table_variantTbl->setRowCount(nrow);
    ui->table_variantTbl->setColumnCount(4);

    for (int i = 0; i < nrow ; i++){

        ui->table_variantTbl->setItem(i, 0, new QTableWidgetItem(QString::number(variants->at(i).trueMaf)));
        ui->table_variantTbl->setItem(i, 1, new QTableWidgetItem(QString::number(variants->at(i).P[0])));
        ui->table_variantTbl->setItem(i, 2, new QTableWidgetItem(QString::number(variants->at(i).P[1])));
        ui->table_variantTbl->setItem(i, 3, new QTableWidgetItem(QString::number(variants->at(i).P[2])));
    }
}


void TableDisplayWindow::on_table_variantTbl_cellClicked(int row, int column){
    fillGenotypeTable(row);
}
