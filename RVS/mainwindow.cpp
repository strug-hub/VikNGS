#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QThread>
#include <QtConcurrent/QtConcurrent>

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

void MainWindow::on_vcfBrowse_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr("VCF File (*.vcf);;All files (*.*)"));

    if(!fileName.isNull())
        ui->vcfDir->setText(fileName);
}

void MainWindow::on_runSimulationButton_clicked()
{
    SimulationRequest request =
            newSimulationRequest(
                ui->nPopulation->text().toStdString(),
                ui->prevelencePopulation->text().toStdString(),
                ui->nSNP->text().toStdString(),
                ui->meanErrorRate->text().toStdString(),
                ui->errorSD->text().toStdString(),
                ui->oddsRatio->text().toStdString(),
                ui->mafLower->text().toStdString(),
                ui->mafUpper->text().toStdString()
                );


    request.print();

    QFuture<void> future = QtConcurrent::run(startSimulation, request);


}

void MainWindow::on_addGroupButton_clicked()
{

    bool isHrg = ui->highRdoButton->isChecked();
    bool isCase = ui->caseRdoButton->isChecked();

    std::string n = ui->nSample->text().toStdString();
    std::string mean = ui->groupMean->text().toStdString();
    std::string sd = ui->groupSD->text().toStdString();

    int index = ui->groupTable->rowCount();

    ui->groupTable->insertRow(index);
    ui->groupTable->setItem(index, 0, new QTableWidgetItem(n));
    ui->groupTable->setItem(index, 1, new QTableWidgetItem(isCase ? "case" : "control"));
    ui->groupTable->setItem(index, 2, new QTableWidgetItem(isHrg ? "high" : "low"));
    ui->groupTable->setItem(index, 3, new QTableWidgetItem(mean));
    ui->groupTable->setItem(index, 4, new QTableWidgetItem(sd));

}
