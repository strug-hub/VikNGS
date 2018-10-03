#include "MainWindow.h"
#include "SimPlotWindow.h"
#include "ui_mainwindow.h"

const QString sep = ":";

SimulationRequest MainWindow::qConstructRequest(std::vector<SimulationRequestGroup> groups){

    SimulationRequest request;
    request.groups = groups;

    request.testStatistic = Statistic::COMMON;
    request.useBootstrap = false;

    if(ui->qsim_testRareSkatBtn->isChecked()){
        request.testStatistic = Statistic::SKAT;
        request.collapse = ui->qsim_collapseTxt->text().toInt();
        request.useBootstrap = true;

    }
    else if(ui->qsim_testRareCastBtn->isChecked()){
        request.testStatistic = Statistic::CAST;
        request.collapse = ui->qsim_collapseTxt->text().toInt();
        request.useBootstrap = true;
    }

    request.stopEarly = false;
    if (request.useBootstrap){
        request.nboot = ui->qsim_testBootTxt->text().toInt();
        request.stopEarly = ui->qsim_testStopChk->isChecked();
    }
    else
        request.nboot=0;

    request.nsnp = ui->qsim_variantSizeTxt->text().toInt();
    request.r2 = ui->qsim_R2Txt->text().toDouble();

    request.mafMin = ui->qsim_variantMafMinTxt->text().toDouble();
    request.mafMax = ui->qsim_variantMafMaxTxt->text().toDouble();

    return request;
}

void MainWindow::on_qsim_stopBtn_clicked(){
    stopJob();
}

void MainWindow::on_qsim_runBtn_clicked(){
    if(!ui->qsim_runBtn->isEnabled())
        return;

    greyOutput();
    disableRun();

    bool ok=false;

    bool multitest = false;
    for(int i = 0; i< ui->qsim_groupTbl->rowCount(); i++)
        if(ui->qsim_groupTbl->item(i, 0)->text().contains(sep))
            multitest = true;

    try{
        int steps = 1;
        if(multitest){
            ok=false;
            steps = ui->qsim_powerStepTxt->text().toInt(&ok);
            if(!ok || steps < 1)
                throwError(ERROR_SOURCE, "Invalid step value. Expected # steps to be an integer greater than 1.", ui->qsim_powerStepTxt->text().toStdString());
        }

        ok=false;
        int highLow = ui->qsim_groupHighLowTxt->text().toInt(&ok);

        if(!ok || highLow <= 0)
            throwError(ERROR_SOURCE, "Expected high-low read depth cutoff to be an integer greater than 0", std::to_string(highLow));

        std::vector<SimulationRequestGroup> groups = constructGroups(steps, ui->qsim_groupTbl, highLow, Family::NORMAL);
        SimulationRequest request = qConstructRequest(groups);
        request.family = Family::NORMAL;

        ok=false;
        double normalMean = ui->qsim_meanTxt->text().toDouble(&ok);
        if(!ok)
            throwError(ERROR_SOURCE, "Quantitative mean should be a numeric value.", ui->qsim_meanTxt->text().toStdString());
        ok=false;
        double normalSd = ui->qsim_sdTxt->text().toDouble(&ok);
        if(!ok)
            throwError(ERROR_SOURCE, "Quantitative SD should be a numeric value.", ui->qsim_sdTxt->text().toStdString());

        for(size_t g = 0; g < request.groups.size(); g++){
            request.groups[g].normalMean = normalMean;
            request.groups[g].normalSd = normalSd;
        }

        request.steps = steps;

        int nthreads = ui->qsim_threadsTxt->text().toInt(&ok);
        if(!ok)
            throwError(ERROR_SOURCE, "Invalid number of threads. Expected an integer greater than or equal to 1.", ui->qsim_threadsTxt->text().toStdString());

        request.nthreads = nthreads;

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

void MainWindow::on_qsim_groupRemoveBtn_clicked()
{
    QList<QPersistentModelIndex> indexes;

    foreach (const QModelIndex &i, ui->qsim_groupTbl->selectionModel()->selectedIndexes())
        indexes << i;

    foreach (const QPersistentModelIndex &i, indexes)
        ui->qsim_groupTbl->model()->removeRow(i.row());


    int index = ui->qsim_groupTbl->rowCount();
    ui->qsim_groupTbl->selectRow(index-1);
}

void MainWindow::qSimEnableRare(bool value){
    ui->qsim_collapseTxt->setEnabled(value);
    ui->qsim_collapseLbl->setEnabled(value);
    ui->qsim_testBootLbl->setEnabled(value);
    ui->qsim_testBootTxt->setEnabled(value);
    ui->qsim_testStopChk->setEnabled(value);
}

void MainWindow::on_qsim_testRareCastBtn_toggled(bool checked){
    qSimEnableRare(checked);
}

void MainWindow::on_qsim_testRareSkatBtn_toggled(bool checked){
    qSimEnableRare(checked);
}

void MainWindow::on_qsim_groupAddBtn_clicked(){
    addGroup(ui->qsim_groupTbl, "100", "normal", "10", "2", "0.01");
}



