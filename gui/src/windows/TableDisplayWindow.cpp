#include "TableDisplayWindow.h"
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

void TableDisplayWindow::initialize(QString title, Data* data, std::vector<int> testIndexes){

    this->testIndexes = testIndexes;

    nsamples = -1;
    for(size_t i = 0; i < testIndexes.size(); i++)
        if(data->tests[testIndexes[i]].getSampleSize() > nsamples)
            this->nsamples = data->tests[testIndexes[i]].getSampleSize();

    this->setWindowTitle("Table - " + title);
    this->variants = &data->variants;

    this->g = data->sampleInfo.getG();
    this->y = data->sampleInfo.getY();
    if(nsamples > 0){
        VectorXi G = g.block(0, 0, nsamples, g.cols()); this->g = G;
        VectorXd Y = y.block(0, 0, nsamples, y.cols()); this->y = Y;
    }
    this->tests = &data->tests;
    buildVariantTable();
}

void TableDisplayWindow::fillGenotypeTable(int variantIndex){

    ui->table_genoTbl->setColumnCount(0);

    int nrow = y.rows();

    Variant* variant;
    int i = 0;
    for(size_t n = 0; n < variants->size(); n++){
        for(size_t m = 0; m < variants->at(n).size(); m++){
            if(i == variantIndex){
                variant = &variants->at(n).getVariants()->at(m);
                goto LOOP_BREAK;
            }
            i++;
        }
    }
    return;
    LOOP_BREAK:

    QStringList titles;
    QVector<QVector<QTableWidgetItem*>> table;

    addSampleInfo(nrow, titles, table);
    addSampleGenotypes(nrow, titles, table, variant);

    int ncol = titles.size();
    ui->table_genoTbl->setColumnCount(ncol);
    ui->table_genoTbl->setRowCount(nrow);
    ui->table_genoTbl->setHorizontalHeaderLabels(titles);

    for(int i = 0; i < table.size(); i++)
        for(int j = 0; j < table[i].size(); j++)
            ui->table_genoTbl->setItem(j,i,table[i][j]);

    ui->table_genoTbl->setHorizontalHeaderLabels(titles);
}

void TableDisplayWindow::buildVariantTable(){

    ui->table_variantTbl->clearContents();
    int nrow = 0;
    for(size_t i = 0; i < variants->size(); i++)
        nrow += variants->at(i).size();

    QStringList titles;
    QVector<QVector<QTableWidgetItem*>> table;

    addInfo(nrow, titles, table);
    addPvals(nrow, titles, table);

    //if caseControl
    addMafsCaseControl(nrow, titles, table, true);
    addMafsCaseControl(nrow, titles, table, false);

    int ncol = titles.size();
    ui->table_variantTbl->setColumnCount(ncol);
    ui->table_variantTbl->setRowCount(nrow);
    ui->table_variantTbl->setHorizontalHeaderLabels(titles);

    for(int i = 0; i < table.size(); i++)
        for(int j = 0; j < table[i].size(); j++)
            ui->table_variantTbl->setItem(j,i,table[i][j]);
}

void TableDisplayWindow::on_table_variantTbl_cellClicked(int row, int column){

    int index = ui->table_variantTbl->item(row,0)->text().toInt();
    fillGenotypeTable(index);
/*
    QString ref = QString::fromStdString(variants->at(index).ref);
    QString alt = QString::fromStdString(variants->at(index).alt);
    ui->table_refLbl->setText(basePair(ref));
    ui->table_altLbl->setText(basePair(alt));
    selectedVariantIndex = index;
*/
}
