#include "mainwindow.h"
#include "plotwindow.h"

#include "ui_mainwindow.h"
#include <iostream>
#include "../log/qlog.h"
#include "../log/typeconverter.h"
#include "../src/Log.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    setWindowIcon(QIcon(":icon.svg"));
    //required to push log updates to textbox
    connect(getQLog(), SIGNAL(pushOutput(QString, QColor)), this, SLOT(printOutput(QString, QColor)));
    qRegisterMetaType<QVector<Variant> >("QVector<Variant>");
    qRegisterMetaType<std::vector<std::vector<Variant>>>("std::vector<std::vector<Variant>>");
    qRegisterMetaType<std::vector<SimulationRequest>>("std::vector<SimulationRequest>");

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::printOutput(QString str, QColor c){

    QString oldOut = ui->outputBox->toHtml();
    QString newOut = "<font color=\"" + c.name() + "\">" + str + "</font>";

    //avoid having the first line blank
    if(ui->outputBox->toPlainText().size() < 2)
        ui->outputBox->setHtml(newOut);
    else
        ui->outputBox->setHtml(oldOut + newOut);

    ui->outputBox->verticalScrollBar()->setValue(ui->outputBox->verticalScrollBar()->maximum());
    ui->outputBox->repaint();
}

void MainWindow::greyOutput(){

    //return if nothing to grey
    if(ui->outputBox->toPlainText().size() < 2)
        return;

    QString oldOut = ui->outputBox->toHtml();
    oldOut = oldOut.replace(QRegExp("color:#.{6};"), "color:" + grey.name() + ";");

    ui->outputBox->setHtml(oldOut + "<br>");
    ui->outputBox->verticalScrollBar()->setValue(ui->outputBox->verticalScrollBar()->maximum());
    ui->outputBox->repaint();

}

void MainWindow::jobFinished(QVector<Variant> variants){
    ui->main_runBtn->setEnabled(true);

    if(variants.size() > 0){
        PlotWindow *plotter = new PlotWindow();
        QString title = "Plot " + QString::number(plotCount);
        plotCount++;
        plotter->initialize(variants, title);
        printOutput("Displaying results in " + title.toLower(), green);
        plotter->show();
    }
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

    greyOutput();

    if(!ui->main_runBtn->isEnabled())
        return;

    ui->main_runBtn->setEnabled(false);

    try{

        QVector<std::string> commands;
        commands.push_back("vikNGS");

        printInfo("Starting vikNGS...");

        std::string vcfDir = ui->main_vcfDirTxt->text().toStdString();
        printInfo("VCF file: " + vcfDir);
        //commands.push_back("--vcf " + vcfDir);
        commands.push_back("--vcf [...]");

        std::string sampleDir = ui->main_sampleDirTxt->text().toStdString();
        printInfo("Sample info file: " + sampleDir);
        //commands.push_back("--sample " + sampleDir);
        commands.push_back("--sample [...]");

        initializeRequest(vcfDir, sampleDir);

        std::string bedDir = ui->main_bedDirTxt->text().toStdString();

        bool collapseGene = ui->main_bedCollapseGeneBtn->isChecked();
        bool collapseCoding = ui->main_bedCollapseCodingBtn->isChecked();
        bool collapseExon = ui->main_bedCollapseCodingBtn->isChecked();

        if(bedDir.size() > 0){

            setCollapseFile(bedDir);
            printInfo("BED file: " + bedDir);
            //commands.push_back("-b " + bedDir);
            commands.push_back("-b [...]");

            if(collapseGene){
                setCollapseGene();
                printInfo("Collapse variants along genes");
            }
            else if(collapseCoding){
                setCollapseCoding();
                printInfo("Collapse variants along coding regions");
            }
            else if(collapseExon){
                setCollapseExon();
                printInfo("Collapse variants along exons");
            }

        }


        std::string mafCutOff = ui->main_vcfMafTxt->text().toStdString();
        double maf = toDouble("Minor allele frequency threshold", mafCutOff);
        setMAFCutoff(maf);
        mafCutOff = toString(maf);
        commands.push_back("-m " + mafCutOff);
        printInfo("Minor allele frequency threshold: " + mafCutOff);

        std::string highLowCutOff = ui->main_sampleDepthTxt->text().toStdString();
        int depth = std::stoi(highLowCutOff);
        setHighLowCutOff(depth);
        highLowCutOff = toString(depth);
        commands.push_back("-d " + highLowCutOff);
        printInfo("Read depth threshold: " + highLowCutOff);

        std::string missingThreshold = ui->main_vcfMissingTxt->text().toStdString();
        double missing = toDouble("Missing data threshold", missingThreshold);
        setMissingThreshold(missing);
        missingThreshold = toString(missing);
        commands.push_back("-x " + missingThreshold);
        printInfo("Missing data threshold: " + missingThreshold);

        bool mustPass = ui->main_vcfPassChk->isChecked();
        setMustPASS(mustPass);
        if(!mustPass){
            commands.push_back("-a");
            printInfo("Retain variants which do not PASS");
        }
        else
            printInfo("Remove variants which do not PASS");

        std::string test = "common";
        if(ui->main_testRareCastBtn->isChecked())
            test = "cast";
        else if(ui->main_testRareCalphaBtn->isChecked())
            test = "calpha";

        if(test=="common"){
          useCommonTest();
          printInfo("Preparing to run common variant association...");
          commands.push_back("-c");
        }
        else if(test=="calpha"){
          useRareTest(test);
          printInfo("Preparing to run rare variant association (C-alpha p-values)...");
          commands.push_back("-r calpha");
        }
        else if(test=="cast"){
          useRareTest(test);
          printInfo("Preparing to run rare variant association (CAST p-values)...");
          commands.push_back("-r cast");
        }

        std::string nboot = ui->main_testBootTxt->text().toStdString();
        bool stopEarly = ui->main_testStopChk->isChecked();

        int n = toInt("Bootstrap iterations", nboot);
        nboot = toString(n);

        if(n > 1 && ui->main_testBootChk->isChecked()){
            useBootstrap(n);

            setStopEarly(stopEarly);
            if(stopEarly){
                printInfo("Using " + nboot + " bootstrap iterations and early stopping");
                commands.push_back("-s");
            }
            else
                printInfo("Using " + nboot + " bootstrap iterations");

            commands.push_back("-n " + nboot);
        }

        std::string nthreads = ui->main_threadsTxt->text().toStdString();
        int threads = toInt("Number of threads", nthreads);
        nthreads = toString(threads);
        setNumberThreads(threads);
        if(threads != 1){
            printInfo("Using " + nthreads + " threads");
            commands.push_back("-t " + nthreads);
        }

        /*

        //TODO:
        //if(k->count() > 0 && b->count() < 1){
        //  setCollapse(collapse);
        //  printInfo("Collapsing every " + std::to_string(collapse) + " variants");
        //}

        //TODO:
        //  if(from != -1 && to != -1){
        //    setMinPos(from);
        //    setMaxPos(to);
        //    printInfo("Analyzing variants between POS " + std::to_string(from) + " and " + std::to_string(to));
        //  }
        //  else if(from != -1){
        //    setMinPos(from);
        //    printInfo("Analyzing variants with POS greater than " + std::to_string(from));
        //  }
        //  else if(to != -1){
        //    setMinPos(to);
        //    printInfo("Analyzing variants with POS less than " + std::to_string(to));
        //  }








*/
    Request req = getRequest();

    QString command = "";
    for (int i = 0; i<commands.size(); i++)
        command.append(QString::fromStdString(commands[i]) + " ");

    printOutput("\n---------------------", green);
    printOutput("Command Line :", green);
    printOutput(command, green);
    printOutput("---------------------\n", green);

    QThread* thread = new QThread;
    Runner* runner = new Runner();
    runner->setRequest(req);
    runner->moveToThread(thread);
    connect(thread, SIGNAL(started()), runner, SLOT(runVikngs()));
    connect(runner, SIGNAL(complete()), thread, SLOT(quit()));
    connect(runner, SIGNAL(complete()), runner, SLOT(deleteLater()));
    connect(runner, SIGNAL(jobFinished(QVector<Variant>)), this, SLOT(jobFinished(QVector<Variant>)));

    thread->start();

    }catch(...){
        ui->main_runBtn->setEnabled(true);
    }
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
    ui->main_testStopChk->setEnabled(ui->main_testBootChk->isChecked());

}

void MainWindow::on_main_testBootChk_toggled(bool checked){

    if(!ui->main_testCommonBtn->isChecked() && !checked){
        ui->main_testBootChk->setChecked(true);
        ui->main_testStopChk->setEnabled(true);
    }

}

