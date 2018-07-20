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

void TableDisplayWindow::initialize(QString title, Result* result){

    simulation = false;

    this->setWindowTitle("Table - " + title);
    this->variants = &result->variants;
    this->g = result->input.G;
    this->y = result->input.Y;

    buildVariantTable();
}

void TableDisplayWindow::drawVariantCheckTree(){

    for(CheckTree* tree : variantCheckTree){
        QTreeWidgetItem* branch = tree->draw();
        ui->table_checkTre->addTopLevelItem(branch);
        connect(ui->table_checkTre, SIGNAL(itemSelectionChanged()),
                tree, SLOT(updateState()));
    }
}

void TableDisplayWindow::fillGenotypeTable(int variantIndex){

    ui->table_genoTbl->setColumnCount(0);

    int nrow = y.rows();

    Variant variant = variants->at(variantIndex);

    QStringList titles;
    QVector<QVector<QTableWidgetItem*>> table;

    addSampleInfo(nrow, titles, table, variant);
    addSampleGenotypes(nrow, titles, table, variant);
    addVCFData(nrow, titles, table, variant);
    addAlleleProbability(nrow, titles, table, variant);

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
    if(simulation && currentRow >= 0){
        int index = ui->table_genoTbl->item(currentRow,0)->text().toInt();
        displayCalls(index);
    }
}
