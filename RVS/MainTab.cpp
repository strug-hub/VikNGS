#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "mainwindow.h"
#include "ui_mainwindow.h"

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


void MainWindow::on_main_vcfDirBtn_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr("VCF File (*.vcf);;All files (*.*)"));
    if(!fileName.isNull())
        ui->main_vcfDirTxt->setText(fileName);
}

void MainWindow::on_main_sampleDirBtn_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr("Text File (*.txt);;All files (*.*)"));
    if(!fileName.isNull())
        ui->main_sampleDirTxt->setText(fileName);
}

void MainWindow::on_main_bedDirBtn_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr("BED File (*.bed);;Text File (*.txt);;All files (*.*)"));
    if(!fileName.isNull())
        ui->main_bedDirTxt->setText(fileName);
}

void MainWindow::on_main_runBtn_clicked()
{
    if(!ui->main_runBtn->isEnabled())
        return;

    std::cout<< "RUN BABY RUN\n\n------\n";

    ui->main_runBtn->setEnabled(false);

    std::string vcfDir = ui->main_vcfDirTxt->text().toStdString();
    std::string sampleDir = ui->main_sampleDirTxt->text().toStdString();
    std::string bedDir = ui->main_bedDirTxt->text().toStdString();

    std::string highLowCutOff = ui->main_sampleDepthTxt->text().toStdString();
    bool collapseCoding = ui->main_bedCollapseCodingBtn->isChecked();
    bool collapseExon = ui->main_bedCollapseCodingBtn->isChecked();

    std::string mafCutOff = ui->main_vcfMafTxt->text().toStdString();
    std::string missingThreshold = ui->main_vcfMissingTxt->text().toStdString();

    bool onlySnps = true;
    bool mustPass = ui->main_vcfPassChk->isChecked();

    std::string test = "common";
    if(ui->main_testRareCastBtn->isChecked())
        test = "cast";
    else if(ui->main_testRareCalphaBtn->isChecked())
        test = "calpha";

    bool boot = ui->main_testBootChk->isChecked();
    std::string nboot = ui->main_testBootTxt->text().toStdString();

    try{
    Request request =
            newRequest(
                vcfDir, sampleDir, bedDir,
                highLowCutOff, collapseCoding, collapseExon,
                mafCutOff, missingThreshold, onlySnps, mustPass,
                test, boot, nboot
                );

    QFuture<std::vector<double>> future = QtConcurrent::run(startVikNGS, request);

    while(future.isRunning()){  }

    QVector<double> pval = QVector<double>::fromStdVector(future.result());

    main_replot(pval);

    }catch(std::invalid_argument e){
        QMessageBox *dialog = new QMessageBox;
        dialog->setWindowTitle("Input error");

        QString error(e.what());


        dialog->setText("Invalid input. Expected " +
                        error.split(",").at(1) +
                        " for following value: " +
                        error.split(",").at(0) + ".");
        dialog->show();

        ui->main_runBtn->setEnabled(true);

        return;

    }catch(std::domain_error e){
        QMessageBox *dialog = new QMessageBox;
        dialog->setWindowTitle("Invalid input");

        QString error(e.what());


        dialog->setText("Invalid input. Please change the following input: \n\n" + error);
        dialog->show();

        ui->sim_runBtn->setEnabled(true);

        return;
    }

    ui->main_runBtn->setEnabled(true);

}


void MainWindow::main_replot(QVector<double> values){

    int n = values.size();

    if(n < 0)
        return;

    QVector<double> x(n), y(n);
    QVector<double> x2(n), y2(n);

    int count = 0;
    std::vector<double> percent;
    percent.push_back(24895.0/308828);
    percent.push_back(24219.0/308828);
    percent.push_back(19829.0/308828);
    percent.push_back(19021.0/308828);
    percent.push_back(18153.0/308828);
    percent.push_back(17080.0/308828);
    percent.push_back(15934.0/308828);
    percent.push_back(14513.0/308828);
    percent.push_back(13839.0/308828);
    percent.push_back(13379.0/308828);
    percent.push_back(13508.0/308828);
    percent.push_back(13327.0/308828);
    percent.push_back(11436.0/308828);
    percent.push_back(10704.0/308828);
    percent.push_back(10199.0/308828);
    percent.push_back(9033.0/308828);
    percent.push_back(8325.0/308828);
    percent.push_back(8037.0/308828);
    percent.push_back(5861.0/308828);
    percent.push_back(6444.0/308828);
    percent.push_back(4670.0/308828);
    percent.push_back(5081.0/308828);
    percent.push_back(15604.0/308828);

    for(int j = 0; j < 23; j++){

        int temp = 0;

        for (int i=count; i < percent[j]*(1.0*n) + count; ++i)
        {
           if(j % 2 == 0){
              x[i] = i;
              y[i] = -log10(values.at(i));
           }
           else{
               x2[i] = i;
               y2[i] = -log10(values.at(i));
           }
           temp = i;
        }
        count = temp;

    }

    ui->main_plotBox->clearPlottables();
    ui->main_plotBox->addGraph();

    // give the axes some labels:
    ui->main_plotBox->xAxis->setLabel("Position");
    ui->main_plotBox->yAxis->setLabel("-log(p)");

    ui->main_plotBox->graph()->setData(x, y);
    ui->main_plotBox->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);
    ui->main_plotBox->graph()->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                Qt::gray, Qt::white, 2));

    ui->main_plotBox->addGraph();

    ui->main_plotBox->graph()->setData(x2, y2);
    ui->main_plotBox->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);
    ui->main_plotBox->graph()->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                Qt::darkGray, Qt::white, 2));

    // set axes ranges, so we see all data:
    ui->main_plotBox->xAxis->setRange(0, n);
    ui->main_plotBox->yAxis->setRange(0, 12);
    ui->main_plotBox->replot();
}

void MainWindow::on_main_testRareCastBtn_toggled(bool checked){
    if(checked)
        ui->main_testBootChk->setChecked(true);
}

void MainWindow::on_main_testRareCalphaBtn_toggled(bool checked){
    if(checked)
        ui->main_testBootChk->setChecked(true);
}

void MainWindow::on_main_testBootChk_stateChanged(int arg1)
{
    ui->main_testBootTxt->setEnabled(ui->main_testBootChk->isChecked());
}

void MainWindow::on_main_testBootChk_toggled(bool checked){

    if(!ui->main_testCommonBtn->isChecked() && !checked)
        ui->main_testBootChk->setChecked(true);

}

