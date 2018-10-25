#include "MainWindow.h"
#include "SimPlotWindow.h"
#include "ui_mainwindow.h"

const QString sep = ":";

void MainWindow::simulationTabInit(){
    ui->sim_simulationSld->setValue(0);
    switchToCaseControl();
    addGroup(ui->sim_groupTbl, "100", "case", "32", "8", "0.01");
    addGroup(ui->sim_groupTbl, "100:300", "control", "6", "2", "0.01");
    prevR2 = "0.0";
    prevOddsRatio = "1.0";

}

void MainWindow::on_sim_runBtn_clicked(){
    if(!ui->sim_runBtn->isEnabled())
        return;

    greyOutput();
    disableRun();

    try{
        SimulationRequest request;

        request.testStatistic = getTestSim();
        int nboot = getBootstrapIterationsSim();
        request.nboot = nboot;
        request.useBootstrap = nboot > 0;
        request.collapse = getCollapseSizeSim();
        request.stopEarly = getEarlyStopSim();

        request.nsnp = ui->sim_variantSizeTxt->text().toInt();
        request.effectSize = ui->sim_effectSizeTxt->text().toDouble();
        request.mafMin = ui->sim_variantMafMinTxt->text().toDouble();
        request.mafMax = ui->sim_variantMafMaxTxt->text().toDouble();

        request.steps = getNumberOfStepsSim();
        request.family = getFamilySim();
        request.nthreads = getNumberOfThreadsSim();
        int highLow = getHighLowCutoffSim();

        request.groups = constructGroups(request.steps, highLow, request.family);

        //throws error if invalid
        request.validate();

        QThread* thread = new QThread;
        AsyncJob* runner = new AsyncJob();

        runner->setSimulationRequest(request);
        runner->moveToThread(thread);
        connect(thread, SIGNAL(started()), runner, SLOT(runSimulation()));
        connect(runner, SIGNAL(complete()), thread, SLOT(quit()));
        connect(runner, SIGNAL(complete()), runner, SLOT(deleteLater()));
        connect(runner, SIGNAL(simulationFinished(Data, SimulationRequest)),
                this, SLOT(simulationFinished(Data, SimulationRequest)));
        thread->start();
        jobThread = thread;
    }catch(...){
        ui->sim_runBtn->setEnabled(true);
        enableRun();
    }
}

void MainWindow::on_sim_stopBtn_clicked(){
    STOP_RUNNING_THREAD = true;
    while(!jobThread->isFinished()){
        jobThread->quit();
    }
    greyOutput();
    printInfo("Job stopped.");
    enableRun();
    STOP_RUNNING_THREAD=false;
}

SimulationRequestGroup readRow(QTableWidget* table, int index, int nsteps){
    SimulationRequestGroup g;

    QString sampleSize = table->item(index, 0)->text();
    if(sampleSize.contains(sep)){

        bool ok1 = false;
        bool ok2 = false;

        int nMin = sampleSize.split(sep).at(0).toInt(&ok1);
        int nMax = sampleSize.split(sep).at(1).toInt(&ok2);
        if(!ok1 || !ok2 || nMin <= 0 || nMax <= 0)
            throwError(ERROR_SOURCE, "Invalid group settings. Ensure number of individuals is an integer greater than 0. If you wish to use a range please separate two integers by a colon \":\"");

        if(nMax < nMin){
            int temp = nMax;
            nMax = nMin;
            nMin = temp;
        }

        g.n=nMin;
        g.n_increment = (1.0*(nMax - nMin))/(1.0*std::max(1, nsteps-1));
    }
    else{
        bool ok = false;
        int n = sampleSize.toInt(&ok);
        if(!ok || n <= 0)
            throwError(ERROR_SOURCE, "Invalid group settings. Ensure number of individuals is an integer greater than 0. If you wish to use a range please separate two integers by a colon \":\"");

        g.n = n;
        g.n_increment=0;
    }

    QString cohort = ((QComboBox*)table->indexWidget(table->model()->index(index, 1)))->currentText();
    g.isCase = (cohort == "case");

    bool ok = false;
    QString meanDepth = table->item(index, 2)->text();
    g.meanDepth = meanDepth.toDouble(&ok);
    if(!ok || g.meanDepth <= 0)
        throwError(ERROR_SOURCE, "Invalid group settings. Expected mean read depth to be a numeric value greater than 0.", std::to_string(g.meanDepth));

    ok = false;
    QString sdDepth = table->item(index, 3)->text();
    g.sdDepth = sdDepth.toDouble(&ok);
    if(!ok || g.sdDepth < 0)
        throwError(ERROR_SOURCE, "Invalid group settings. Expected read depth SD to be a numeric value greater than or equal to 0.", std::to_string(g.sdDepth));

    ok = false;
    QString errorRate = table->item(index, 4)->text();
    g.errorRate = errorRate.toDouble(&ok);
    if(!ok || g.errorRate < 0 || g.errorRate > 1)
        throwError(ERROR_SOURCE, "Invalid group settings. Expected error rate to be a numeric value between 0 and 1", std::to_string(g.errorRate));

    return g;
}

std::vector<SimulationRequestGroup> MainWindow::constructGroups(int nsteps, int highLow, Family family){

    std::vector<SimulationRequestGroup> groups;

    for(int i = 0; i < ui->sim_groupTbl->rowCount(); i++){

        SimulationRequestGroup g = readRow(ui->sim_groupTbl, i, nsteps);
        g.index = i;
        g.family = family;
        (g.meanDepth >= highLow) ? g.readDepth = Depth::HIGH : g.readDepth = Depth::LOW;

        if(family == Family::NORMAL){
            bool ok = false;
            QString phenoMean = ui->sim_meanTxt->text();
            QString phenoSd = ui->sim_sdTxt->text();

            g.normalMean = phenoMean.toDouble(&ok);
            if(!ok)
                throwError(ERROR_SOURCE, "Invalid quantitative phenotype mean value. Expected a numeric value", phenoMean.toStdString());

            ok = false;
            g.normalSd = phenoSd.toDouble(&ok);
            if(!ok)
                throwError(ERROR_SOURCE, "Invalid quantitative phenotype standard deviation. Expected a numeric value", phenoSd.toStdString());
        }


        groups.push_back(g);
    }

    return groups;
}

void MainWindow::simulationFinished(Data results, SimulationRequest req){
    enableRun();

    if(results.size() > 0){

        bool showWindow = false;
        for(VariantSet vs : results.variants){
            for(int i = 0; i < vs.nPvals(); i++){
                if(!std::isnan(vs.getPval(i))){
                    showWindow = true;
                    goto done;
                }
            }
        }

        done:

        if(!showWindow){
            printInfo("Simulation stopped");
            return;
        }
        SimPlotWindow *plotter = new SimPlotWindow();
        QString title = "Plot " + QString::number(plotCount);
        plotCount++;
        plotter->initialize(results, req, title);
        printOutput("Displaying results in " + title.toLower(), green);
        plotter->show();
    }
}

void MainWindow::on_sim_groupRemoveBtn_clicked()
{
    QList<QPersistentModelIndex> indexes;

    foreach (const QModelIndex &i, ui->sim_groupTbl->selectionModel()->selectedIndexes())
        indexes << i;

    foreach (const QPersistentModelIndex &i, indexes)
        ui->sim_groupTbl->model()->removeRow(i.row());


    int index = ui->sim_groupTbl->rowCount();
    ui->sim_groupTbl->selectRow(index-1);

}

void MainWindow::simEnableRare(bool value){
    ui->sim_collapseTxt->setEnabled(value);
    ui->sim_collapseLbl->setEnabled(value);
    ui->sim_testBootLbl->setEnabled(value);
    ui->sim_testBootTxt->setEnabled(value);
    ui->sim_testStopChk->setEnabled(value);
}

void MainWindow::on_sim_testRareCastBtn_toggled(bool checked){
    simEnableRare(checked);
}

void MainWindow::on_sim_testRareSkatBtn_toggled(bool checked){
    simEnableRare(checked);
}

void MainWindow::addGroup(QTableWidget* table, QString n, QString cohort, QString depth, QString sdDepth, QString errorRate){

    table->setSelectionBehavior(QAbstractItemView::SelectRows);

    int index = table->rowCount();

    table->insertRow(index);

    QTableWidgetItem* item;

    item = new QTableWidgetItem(n);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsEditable);
    table->setItem(index, 0, item);

    QComboBox *cohortBox = new QComboBox();
    if(cohort == "control" || cohort == "case"){
        cohortBox->addItem("case");
        cohortBox->addItem("control");

        if(cohort == "control")
            cohortBox->setCurrentIndex(1);
    }
    else{
        cohortBox->addItem("normal");
    }

    table->setIndexWidget(table->model()->index(index, 1), cohortBox);

    item = new QTableWidgetItem(depth);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsEditable);
    table->setItem(index, 2, item);

    item = new QTableWidgetItem(sdDepth);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsEditable);
    table->setItem(index, 3, item);

    item = new QTableWidgetItem(errorRate);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsEditable);
    table->setItem(index, 4, item);

    table->selectRow(index);
    table->update();
}

void MainWindow::on_sim_groupAddBtn_clicked(){
    addGroup(ui->sim_groupTbl, "100", "control", "10", "2", "0.01");
}

void MainWindow::switchToCaseControl(){

    prevR2 = ui->sim_effectSizeTxt->text();
    ui->sim_effectSizeTxt->setText(prevOddsRatio);
    ui->sim_effectSizeLbl->setText("Odds Ratio");

    for(int i =0; i < ui->sim_groupTbl->rowCount(); i++){
        QString newValue = "case";
        if(prevCaseControlValues.size() > i)
            newValue = prevCaseControlValues[i];

        QComboBox *cohortBox = new QComboBox();
        cohortBox->addItem("case");
        cohortBox->addItem("control");

        if(newValue == "control")
            cohortBox->setCurrentIndex(1);

        ui->sim_groupTbl->setIndexWidget(ui->sim_groupTbl->model()->index(i, 1), cohortBox);
        ui->sim_quantGrp->setVisible(false);
    }
}

void MainWindow::switchToQuantitative(){

    prevOddsRatio = ui->sim_effectSizeTxt->text();
    ui->sim_effectSizeTxt->setText(prevR2);
    ui->sim_effectSizeLbl->setText("R²");

    QVector<QString> prevValues;
    for(int i =0; i < ui->sim_groupTbl->rowCount(); i++){
        prevValues.push_back(((QComboBox*) ui->sim_groupTbl->indexWidget(
                                   ui->sim_groupTbl->model()->index(i, 1)))->currentText());
        QComboBox *cohortBox = new QComboBox();
        cohortBox->addItem("normal");
        ui->sim_groupTbl->setIndexWidget(ui->sim_groupTbl->model()->index(i, 1), cohortBox);

    }
    prevCaseControlValues = prevValues;
    ui->sim_quantGrp->setVisible(true);
}

void MainWindow::on_sim_simulationSld_valueChanged(int value){
    if(value == 0)
        switchToCaseControl();
    else if (value == 1)
        switchToQuantitative();
}

void MainWindow::on_sim_caseControlBtn_pressed()
{
    ui->sim_simulationSld->setValue(0);
}

void MainWindow::on_sim_quantitativeBtn_pressed()
{
    ui->sim_simulationSld->setValue(1);
}


int MainWindow::getNumberOfStepsSim(){
    bool multitest = false;
    for(int i = 0; i< ui->sim_groupTbl->rowCount(); i++)
        if(ui->sim_groupTbl->item(i, 0)->text().contains(sep))
            multitest = true;

    if(!multitest)
        return 1;

    bool ok = false;
    int steps = ui->sim_powerStepTxt->text().toInt(&ok);
    if(!ok || steps < 1)
        throwError(ERROR_SOURCE, "Invalid step value. Expected # steps to be an integer greater than 1.", ui->sim_powerStepTxt->text().toStdString());

    return steps;
}

int MainWindow::getHighLowCutoffSim(){
    bool ok = false;
    int highLow = ui->sim_groupHighLowTxt->text().toInt(&ok);

    if(!ok || highLow <= 0)
        throwError(ERROR_SOURCE, "Expected high-low read depth cutoff to be an integer greater than 0", std::to_string(highLow));
    return highLow;
}

Family MainWindow::getFamilySim(){
    if(ui->sim_simulationSld->value() == 0)
        return Family::BINOMIAL;
    else
        return Family::NORMAL;

    return Family::NONE;
}

int MainWindow::getNumberOfThreadsSim(){
    bool ok = false;
    int nthreads = ui->sim_threadsTxt->text().toInt(&ok);
    if(!ok)
        throwError(ERROR_SOURCE, "Invalid number of threads. Expected an integer greater than or equal to 1.", ui->sim_threadsTxt->text().toStdString());

    return nthreads;
}

Statistic MainWindow::getTestSim(){

    if(ui->sim_testRareSkatBtn->isChecked())
        return Statistic::SKAT;
    else if(ui->sim_testRareCastBtn->isChecked())
        return Statistic::CAST;
    else if(ui->sim_testCommonBtn->isChecked())
        return Statistic::COMMON;

    return Statistic::NONE;
}

int MainWindow::getBootstrapIterationsSim(){

    if(!ui->sim_testBootTxt->isEnabled())
        return 0;

    bool ok = false;
    int nboot = ui->sim_testBootTxt->text().toInt(&ok);

    if(!ok || nboot < 1)
        throwError(ERROR_SOURCE, "Invalid number of test iterations. Expected an integer greater than or equal to 1.",  ui->sim_testBootTxt->text().toStdString());

    return nboot;
}

int MainWindow::getCollapseSizeSim(){

    if(!ui->sim_collapseTxt->isEnabled())
        return 1;

    bool ok = false;
    int collapse = ui->sim_collapseTxt->text().toInt(&ok);

    if(!ok || collapse < 1)
        throwError(ERROR_SOURCE, "Invalid collapse size. Expected an integer greater than or equal to 1.",  ui->sim_collapseTxt->text().toStdString());

    return collapse;
}

bool MainWindow::getEarlyStopSim(){

    if(!ui->sim_testStopChk->isEnabled())
        return false;

    return ui->sim_testStopChk->isChecked();
}
