#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QThread>
#include <QtConcurrent/QtConcurrent>
#include <QMessageBox>
#include <QVector>

QVector<double> pval1;
QVector<double> pval2;

QVector<double> pvalMain1;
QVector<double> pvalMain2;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_runSimulationButton_clicked()
{
    if(!ui->runSimulationButton->isEnabled())
        return;

    std::cout<< "RUN BABY RUN\n\n------\n";

    ui->runSimulationButton->setEnabled(false);

    std::vector<SimulationRequestGroup> groups;

    try{
    for(int i = 0; i< ui->groupTable->rowCount(); i++){
        groups.push_back(
                    newSimulationRequestGroup( i,
                        ui->groupTable->item(i, 0)->text().toStdString(),
                        ui->groupTable->item(i, 1)->text().toStdString(),
                        ui->groupTable->item(i, 2)->text().toStdString(),
                        ui->groupTable->item(i, 3)->text().toStdString(),
                        ui->groupTable->item(i, 4)->text().toStdString()
                    )
                );

        groups[i].print();
    }

    SimulationRequest request =
            newSimulationRequest(
                ui->nPopulation->text().toStdString(),
                ui->prevelencePopulation->text().toStdString(),
                ui->nSNP->text().toStdString(),
                ui->meanErrorRate->text().toStdString(),
                ui->errorSD->text().toStdString(),
                ui->oddsRatio->text().toStdString(),
                ui->mafLower->text().toStdString(),
                ui->mafUpper->text().toStdString(),
                groups,
                ui->commonRdoButton->isChecked(),
                ui->bootstrapChkBox->isChecked(),
                ui->nBoot->text().toStdString()
                );

    request.print();

    QFuture<std::vector<std::vector<double>>> future = QtConcurrent::run(startSimulation, request);

    while(future.isRunning()){  }

    auto results = future.result();

    QVector<double> p1;
    QVector<double> p2;

    for(int i =0; i< results.size(); i++){
        p1.push_back(results[i][0]);
        if(results[i].size() > 1)
            p2.push_back(results[i][1]);
    }

    pval1 = p1;
    pval2 = p2;


    replot();

    }catch(std::invalid_argument e){
        QMessageBox *dialog = new QMessageBox;
        dialog->setWindowTitle("Input error");

        QString error(e.what());


        dialog->setText("Invalid input. Expected " +
                        error.split(",").at(1) +
                        " for following value: " +
                        error.split(",").at(0) + ".");
        dialog->show();

        ui->runSimulationButton->setEnabled(true);

        return;
    }catch(std::domain_error e){
        QMessageBox *dialog = new QMessageBox;
        dialog->setWindowTitle("Invalid input");

        QString error(e.what());


        dialog->setText("Invalid input. Please change the following input: \n\n" + error);
        dialog->show();

        ui->runSimulationButton->setEnabled(true);

        return;
    }

    ui->runSimulationButton->setEnabled(true);

}

void MainWindow::on_addGroupButton_clicked()
{
    ui->groupTable->setSelectionBehavior(QAbstractItemView::SelectRows);

    bool isHrg = ui->highRdoButton->isChecked();
    bool isCase = ui->caseRdoButton->isChecked();

    bool isQuantative = ui->quantRdoButton->isChecked();
    QString quantMean = ui->quantMean->text();
    QString quantSD = ui->quantSD->text();

    QString n = ui->nSample->text();
    QString mean = ui->groupMean->text();
    QString sd = ui->groupSD->text();

    int index = ui->groupTable->rowCount();

    ui->groupTable->insertRow(index);
    ui->groupTable->setItem(index, 0, new QTableWidgetItem(n));

    if(!isQuantative)
        ui->groupTable->setItem(index, 1, new QTableWidgetItem(isCase ? "case" : "control"));
    else{
        ui->groupTable->setItem(index, 1, new QTableWidgetItem(quantMean + " : " + quantSD));
    }

    ui->groupTable->setItem(index, 2, new QTableWidgetItem(isHrg ? "high" : "low"));
    ui->groupTable->setItem(index, 3, new QTableWidgetItem(mean));
    ui->groupTable->setItem(index, 4, new QTableWidgetItem(sd));


    ui->groupTable->selectRow(index);

}

void MainWindow::on_removeGroupButton_clicked()
{
    QList<QPersistentModelIndex> indexes;

    foreach (const QModelIndex &i, ui->groupTable->selectionModel()->selectedIndexes())
        indexes << i;

    foreach (const QPersistentModelIndex &i, indexes)
        ui->groupTable->model()->removeRow(i.row());


    int index = ui->groupTable->rowCount();
    ui->groupTable->selectRow(index-1);

}

void MainWindow::replot(){

    int binFactor = ui->binSlider->value();

    if(pval1.size() < 0)
        return;

    int nbins = binFactor*2;
    QVector<double> x(nbins),
            y1(nbins), y2(nbins);

    for (int i=0; i<nbins; ++i)
    {
      x[i] = (i*1.0)/nbins + 0.5/nbins;
      y1[i] = 0;
      y2[i] = 0;
    }

    for (int i=0; i<pval1.size(); ++i)
    {
        int binIndex = floor(pval1[i]/(1.0/nbins));
        y1[binIndex]++;
    }

    for (int i=0; i<pval2.size(); ++i)
    {
        int binIndex = floor(pval2[i]/(1.0/nbins));
        y2[binIndex]++;
    }

    ui->plotBox->clearPlottables();

    // create graph and assign data to it:

    // give the axes some labels:
    ui->plotBox->xAxis->setLabel("x");
    ui->plotBox->yAxis->setLabel("y");

    QCPBars *bars1 = new QCPBars(ui->plotBox->xAxis, ui->plotBox->yAxis);
    bars1->setData(x, y1);
    bars1->setWidth(1.0/nbins);

    QCPBars *bars2 = new QCPBars(ui->plotBox->xAxis, ui->plotBox->yAxis);
    bars2->setData(x, y2);
    bars2->setWidth(1.0/nbins);
    bars2->setBrush(QColor(250, 170, 20, 75).lighter(150));
    bars2->setPen(QColor(250, 170, 20));

    // set axes ranges, so we see all data:
    ui->plotBox->xAxis->setRange(0, 1.05);
    ui->plotBox->yAxis->setRange(0, pval1.size());
    ui->plotBox->replot();
}

void MainWindow::on_binSlider_valueChanged(int value)
{
    replot();
}


void MainWindow::on_rareRdoButton_clicked()
{
    ui->bootstrapChkBox->setChecked(true);
}


void MainWindow::on_bootstrapChkBox_stateChanged(int arg1)
{
    ui->nBoot->setEnabled(ui->bootstrapChkBox->isChecked());
}


//****************************************************************************************************************************************************************

void MainWindow::on_vcfBrowse_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr("VCF File (*.vcf);;All files (*.*)"));
    if(!fileName.isNull())
        ui->vcfDir->setText(fileName);
}

void MainWindow::on_sampleBrowse_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr("Text File (*.txt);;All files (*.*)"));
    if(!fileName.isNull())
        ui->sampleDir->setText(fileName);
}

void MainWindow::on_sampleBrowse_2_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr("BED File (*.bed);;Text File (*.txt);;All files (*.*)"));
    if(!fileName.isNull())
        ui->bedDir->setText(fileName);
}

void MainWindow::on_runButton_clicked()
{
    if(!ui->runSimulationButton->isEnabled())
        return;

    std::cout<< "RUN BABY RUN\n\n------\n";

    ui->runButton->setEnabled(false);

    try{
    Request request =
            newRequest(
                ui->vcfDir->text().toStdString(),
                ui->sampleDir->text().toStdString(),
                ui->bedDir->text().toStdString(),
                ui->highLowCutOff->text().toStdString(),
                ui->codingCollapseRdoButton->isChecked(),
                ui->exonCollapseRdoButton->isChecked(),
                ui->mafCutOff->text().toStdString(),
                ui->missingThreshold->text().toStdString(),
                true,
                ui->mustPassChkBox->isChecked(),
                ui->mainCommonRdoButton->isChecked(),
                ui->mainBootstrapChkBox->isChecked(),
                ui->mainNBoot->text().toStdString()
                );

    QFuture<std::vector<std::vector<double>>> future = QtConcurrent::run(startVikNGS, request);

    while(future.isRunning()){  }

    auto results = future.result();

    QVector<double> p1;
    QVector<double> p2;

    for(int i =0; i< results.size(); i++){
        p1.push_back(results[i][0]);
        if(results[i].size() > 1)
            p2.push_back(results[i][1]);
    }

    pvalMain1 = p1;
    pvalMain2 = p2;

    replotMain();

    }catch(std::invalid_argument e){
        QMessageBox *dialog = new QMessageBox;
        dialog->setWindowTitle("Input error");

        QString error(e.what());


        dialog->setText("Invalid input. Expected " +
                        error.split(",").at(1) +
                        " for following value: " +
                        error.split(",").at(0) + ".");
        dialog->show();

        ui->runButton->setEnabled(true);

        return;
    }catch(std::domain_error e){
        QMessageBox *dialog = new QMessageBox;
        dialog->setWindowTitle("Invalid input");

        QString error(e.what());


        dialog->setText("Invalid input. Please change the following input: \n\n" + error);
        dialog->show();

        ui->runButton->setEnabled(true);

        return;
    }

    ui->runButton->setEnabled(true);

}

void MainWindow::replotMain(){

    int binFactor = ui->mainBinSlider->value();

    if(pvalMain1.size() < 0)
        return;

    int nbins = binFactor*2;
    QVector<double> x(nbins),
            y1(nbins), y2(nbins);

    for (int i=0; i<nbins; ++i)
    {
      x[i] = (i*1.0)/nbins + 0.5/nbins;
      y1[i] = 0;
      y2[i] = 0;
    }

    for (int i=0; i<pvalMain1.size(); ++i)
    {
        int binIndex = floor(pvalMain1[i]/(1.0/nbins));
        y1[abs(binIndex)]++;
    }

    for (int i=0; i<pvalMain2.size(); ++i)
    {
        int binIndex = floor(pvalMain2[i]/(1.0/nbins));
        y2[binIndex]++;
    }

    std::cout << pvalMain1.size();

    ui->mainPlotBox->clearPlottables();

    // create graph and assign data to it:

    // give the axes some labels:
    ui->mainPlotBox->xAxis->setLabel("x");
    ui->mainPlotBox->yAxis->setLabel("y");

    QCPBars *bars1 = new QCPBars(ui->mainPlotBox->xAxis, ui->mainPlotBox->yAxis);
    bars1->setData(x, y1);
    bars1->setWidth(1.0/nbins);

    QCPBars *bars2 = new QCPBars(ui->mainPlotBox->xAxis, ui->mainPlotBox->yAxis);
    bars2->setData(x, y2);
    bars2->setWidth(1.0/nbins);
    bars2->setBrush(QColor(250, 170, 20, 75).lighter(150));
    bars2->setPen(QColor(250, 170, 20));

    // set axes ranges, so we see all data:
    ui->mainPlotBox->xAxis->setRange(0, 1.05);
    ui->mainPlotBox->yAxis->setRange(0, pval1.size());
    ui->mainPlotBox->replot();
}

void MainWindow::on_mainBinSlider_valueChanged(int value)
{
    replotMain();
}

