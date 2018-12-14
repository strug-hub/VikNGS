#include "MainWindow.h"
#include "PlotWindow.h"

#include "ui_mainwindow.h"
#include <iostream>
#include <random>
#include "../log/qlog.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowIcon(QIcon(":icon.svg"));

    //required to push log updates to textbox
    connect(getQLog(), SIGNAL(pushOutput(QString, QColor)), this, SLOT(printOutput(QString, QColor)));
    qRegisterMetaType<Data>("Data");
    qRegisterMetaType<SimulationRequest>("SimulationRequest");

    simulationTabInit();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::printOutput(QString str, QColor c){

    QString oldOut = ui->outputBox->toHtml();
    oldOut = oldOut.right(50000);
    int close = oldOut.indexOf("<");
    oldOut = oldOut.right(50000 - close);

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

void MainWindow::stopJob(){
    STOP_RUNNING_THREAD = true;
    while(!jobThread->isFinished()){
        jobThread->quit();
    }
    greyOutput();
    printInfo("Job stopped.");
    enableRun();
    STOP_RUNNING_THREAD=false;
}

void MainWindow::enableRun(){
    ui->main_runBtn->setEnabled(true);
    ui->main_stopBtn->setEnabled(false);
    ui->sim_runBtn->setEnabled(true);
    ui->sim_stopBtn->setEnabled(false);
}
void MainWindow::disableRun(){
    ui->main_runBtn->setEnabled(false);
    ui->main_stopBtn->setEnabled(true);
    ui->sim_runBtn->setEnabled(false);
    ui->sim_stopBtn->setEnabled(true);
}

void MainWindow::jobFinished(Data result){
    enableRun();
    if(result.size() > 0){
        PlotWindow *plotter = new PlotWindow();
        QString title = "Plot " + QString::number(plotCount);
        plotCount++;
        plotter->initialize(result, title);
        printOutput("Displaying results in " + title.toLower(), green);
        plotter->show();
    }
}

void MainWindow::on_main_stopBtn_clicked()
{
    STOP_RUNNING_THREAD = true;
    while(!jobThread->isFinished()){
        jobThread->quit();
    }
    greyOutput();
    printInfo("Job stopped.");
    enableRun();
    STOP_RUNNING_THREAD=false;
}

void MainWindow::on_main_bedCollapseKBtn_toggled(bool checked)
{
    if(checked)
        ui->main_bedCollapseKTxt->setText("5");
}

void MainWindow::on_main_bedCollapseExonBtn_toggled(bool checked)
{
    if(checked)
        ui->main_bedCollapseKTxt->setText("100");
}

void MainWindow::on_main_bedCollapseGeneBtn_toggled(bool checked)
{
    if(checked)
        ui->main_bedCollapseKTxt->setText("100");
}


void MainWindow::on_main_randomBtn_pressed(){
    ui->main_randomBtn->setVisible(false);

    PlotWindow *plotter = new PlotWindow();
    QString title = "Random";
    Data result; std::vector<VariantSet> vss;

     for(int h = 1; h <= 23; h++){

         double maxSigFactor = 1;
         double oldSigFactor = 1;
         double currentSigFactor = 1;
         double decayRate = 10;

         int size = 10000;
         double mx = 1e6;
         int maxpos = mx + mx*(22-h)/5 + (mx)* randomDouble(0,1)/5;
         std::vector<int> pos(size);
         for(int i =0; i< size; i++)
             pos[i] = randomInt(1,maxpos);

        std::sort(begin(pos), end(pos), std::greater<int>());

         int prevPos = -1;

         for(int i = 0; i < pos.size(); i++){

             if (prevPos == pos[i])
                     continue;

             if(maxSigFactor < 1.01 && randomDouble(0,1) < 1.0/(pos.size()*5))
                 maxSigFactor = randomInt(4,20);

             if(maxSigFactor > 1){
                 currentSigFactor += maxSigFactor * randomDouble(0,1)/decayRate;
                 if(currentSigFactor >= maxSigFactor){
                     currentSigFactor = maxSigFactor;
                     oldSigFactor = maxSigFactor;
                     maxSigFactor = 1;
                 }
             }

             currentSigFactor -= oldSigFactor * randomDouble(0,1)/decayRate;
             if(currentSigFactor <= 1){
                 currentSigFactor = 1;
                 oldSigFactor = 1;
             }

             VariantSet vs;
             if(randomDouble(0,1) < 0.5)
                vs.addPval(randomDouble(0,1)/pow(10, currentSigFactor-1));
             else
                 vs.addPval(randomDouble(0,1));

             Variant v(std::to_string(h), pos[i], "uid", "A", "T");
             vs.addVariant(v);
             vss.push_back(vs);
             prevPos = pos[i];
         }
     }

     result.variants = vss;
     plotter->initialize(result, title);
     plotter->show();
}

void MainWindow::on_pushButton_pressed()
{
    QDesktopServices::openUrl(QUrl("https://vikngsdocs.readthedocs.io/en/latest/index.html"));

}
