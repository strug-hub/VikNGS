#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QThread>
#include <QtConcurrent/QtConcurrent>
#include <QMessageBox>

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
                groups
                );

    request.print();

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

    //QFuture<void> future = QtConcurrent::run(startSimulation, request);


}



void MainWindow::on_addGroupButton_clicked()
{
    ui->groupTable->setSelectionBehavior(QAbstractItemView::SelectRows);

    bool isHrg = ui->highRdoButton->isChecked();
    bool isCase = ui->caseRdoButton->isChecked();


    QString n = ui->nSample->text();
    QString mean = ui->groupMean->text();
    QString sd = ui->groupSD->text();

    int index = ui->groupTable->rowCount();

    ui->groupTable->insertRow(index);
    ui->groupTable->setItem(index, 0, new QTableWidgetItem(n));
    ui->groupTable->setItem(index, 1, new QTableWidgetItem(isCase ? "case" : "control"));
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

