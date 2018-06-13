#include "mainwindow.h"
#include "simplotwindow.h"

#include "ui_mainwindow.h"

const QString sep = ":";

std::vector<SimulationRequestGroup> MainWindow::constructGroups(int test, int ntest){

    std::vector<SimulationRequestGroup> groups;

    for(int i = 0; i< ui->sim_groupTbl->rowCount(); i++){

        SimulationRequestGroup g;

        QString size = ui->sim_groupTbl->item(i, 0)->text();
        int n;

        if(size.contains(sep)){
            int nMin = size.split(sep).at(0).toInt();
            int nMax = size.split(sep).at(1).toInt();

            double denom = std::max(1, ntest-1);
            n = (nMax - nMin)/denom * test + nMin;
        }
        else
            n = size.toInt();

        QString caseControl = ui->sim_groupTbl->item(i, 1)->text();
        QString meanDepth = ui->sim_groupTbl->item(i, 2)->text();
        QString sdDepth = ui->sim_groupTbl->item(i, 3)->text();
        QString errorRate = ui->sim_groupTbl->item(i, 4)->text();

        QString highLow =
                ui->sim_groupTbl->item(i, 2)->text().toDouble() < ui->sim_groupHighLowTxt->text().toDouble() ?
                    "low" : "high";

        g.index = i;
        g.n = n;
        g.isCase = (caseControl == "case");
        g.isHrg = (highLow == "high");
        g.meanDepth = meanDepth.toDouble();
        g.sdDepth = sdDepth.toDouble();
        g.errorRate = errorRate.toDouble();

        groups.push_back(g);

    }

    return groups;
}

SimulationRequest MainWindow::constructRequest(std::vector<SimulationRequestGroup> groups){

    SimulationRequest request;
    request.groups = groups;

    request.test = "common";
    if(ui->sim_testRareCalphaBtn->isChecked()){
        request.test = "calpha";
        request.collapse = ui->sim_collapseTxt->text().toInt();
    }
    else if(ui->sim_testRareCastBtn->isChecked()){
        request.test= "cast";
        request.collapse = ui->sim_collapseTxt->text().toInt();
    }

    if(ui->sim_rvsChk->isChecked())
        request.rvs = true;
    else
        request.rvs = false;

    if(ui->sim_gtChk->isChecked())
        request.regular = true;
    else
        request.regular = false;

    request.useBootstrap = ui->sim_testBootChk->isChecked();
    request.stopEarly = false;
    if (request.useBootstrap){
        request.nboot = ui->sim_testBootTxt->text().toInt();
        request.stopEarly = ui->sim_testStopChk->isChecked();
    }
    else
        request.nboot=0;

    request.nsnp = ui->sim_variantSizeTxt->text().toInt();
    request.oddsRatio = ui->sim_oddsRatioTxt->text().toDouble();
    request.mafMin = ui->sim_variantMafMinTxt->text().toDouble();
    request.mafMax = ui->sim_variantMafMaxTxt->text().toDouble();

    return request;
}

void MainWindow::simulationFinished(std::vector<std::vector<Variant>> variants, std::vector<SimulationRequest> reqs){
    ui->sim_runBtn->setEnabled(true);

    if(variants.size() > 0){
        SimPlotWindow *plotter = new SimPlotWindow();
        QString title = "Plot " + QString::number(plotCount);
        plotCount++;
        plotter->initialize(variants, reqs, title);
        printOutput("Displaying results in " + title.toLower(), green);
        plotter->show();
    }
}

void MainWindow::on_sim_runBtn_clicked(){

    if(!ui->sim_runBtn->isEnabled())
        return;

    ui->sim_runBtn->setEnabled(false);
    greyOutput();

    int ntest = 1;
    bool multitest = false;
    for(int i = 0; i< ui->sim_groupTbl->rowCount(); i++)
        if(ui->sim_groupTbl->item(i, 0)->text().contains(sep))
            multitest = true;

    if(multitest)
        ntest = std::max(1, ui->sim_powerStepTxt->text().toInt());

    std::vector<SimulationRequest> requests;

    try{

        for(int run = 0; run < ntest; run++){

            std::vector<SimulationRequestGroup> groups = constructGroups(run, ntest);
            SimulationRequest request = constructRequest(groups);

            //throws error if invalid
            request.validate();

            requests.push_back(request);
        }

        QThread* thread = new QThread;
        Runner* runner = new Runner();
        runner->setSimulationRequests(requests);
        runner->moveToThread(thread);
        connect(thread, SIGNAL(started()), runner, SLOT(runSimulation()));
        connect(runner, SIGNAL(complete()), thread, SLOT(quit()));
        connect(runner, SIGNAL(complete()), runner, SLOT(deleteLater()));
        connect(runner, SIGNAL(simulationFinished(std::vector<std::vector<Variant>>, std::vector<SimulationRequest>)),
                this, SLOT(simulationFinished(std::vector<std::vector<Variant>>, std::vector<SimulationRequest>)));

        thread->start();
    }catch(...){
        ui->sim_runBtn->setEnabled(true);
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

void MainWindow::on_sim_testRareCastBtn_toggled(bool checked){
    if(checked)
        ui->sim_testBootChk->setChecked(true);

    if(checked)
        ui->sim_collapseTxt->setEnabled(true);
    else
        ui->sim_collapseTxt->setEnabled(false);

}

void MainWindow::on_sim_testRareCalphaBtn_toggled(bool checked){
    if(checked)
        ui->sim_testBootChk->setChecked(true);

    if(checked)
        ui->sim_collapseTxt->setEnabled(true);
    else
        ui->sim_collapseTxt->setEnabled(false);
}

void MainWindow::on_sim_testBootChk_stateChanged(int arg1)
{
    ui->sim_testBootTxt->setEnabled(ui->sim_testBootChk->isChecked());

}

void MainWindow::on_sim_testBootChk_toggled(bool checked){

    if(!ui->sim_testCommonBtn->isChecked() && !checked)
        ui->sim_testBootChk->setChecked(true);
}


void MainWindow::on_sim_groupAddBtn_clicked()
{
    ui->sim_groupTbl->setSelectionBehavior(QAbstractItemView::SelectRows);

    bool isCase = ui->sim_groupCaseBtn->isChecked();

    QString n = "";

    QString nMin = ui->sim_groupSizeMinTxt->text();
    QString nMax = ui->sim_groupSizeMaxTxt->text();

    bool ok1 = false;
    bool ok2 = false;

    int mn = nMin.toInt(&ok1);
    int mx = nMax.toInt(&ok2);

    if(!ok1 || !ok2 || mn <= 0 || mx <= 0){
        warningDialog("Invalid group settings. Ensure number of individuals is an integer greater than 0.");
        return;
    }

    if(mn > mx){
        QString temp = nMin;
        nMin = nMax;
        nMax = temp;
    }

    n = (nMin == nMax) ? nMin : nMin + sep + nMax;


    QString mean = ui->sim_groupDepthMeanTxt->text();
    QString sd = ui->sim_groupDepthSdTxt->text();
    QString error = ui->sim_groupErrorTxt->text();

    bool ok = false;
    double checkMean = mean.toDouble(&ok);

    if( !ok || checkMean <= 0){
        warningDialog("Invalid group settings. Expected mean read depth to be a numeric value greater than 0.");
        return;
    }

    ok=false;
    double checkSd = sd.toDouble(&ok);

    if(!ok || checkSd < 0){
        warningDialog("Invalid group settings. Expected read depth SD to be a numeric value greater than or equal to 0.");
        return;
    }

    ok=false;
    double checkError = error.toDouble(&ok);

    if(!ok || error < 0){
        warningDialog("Invalid group settings. Expected error rate to be a numeric value greater than or equal to 0.");
        return;
    }

    int index = ui->sim_groupTbl->rowCount();

    ui->sim_groupTbl->insertRow(index);

    QTableWidgetItem* item;

    item = new QTableWidgetItem(n);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);
    ui->sim_groupTbl->setItem(index, 0, item);

    item = new QTableWidgetItem(isCase ? "case" : "control");
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);
    ui->sim_groupTbl->setItem(index, 1, item);

    item = new QTableWidgetItem(mean);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);
    ui->sim_groupTbl->setItem(index, 2, item);

    item = new QTableWidgetItem(sd);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);
    ui->sim_groupTbl->setItem(index, 3, item);

    item = new QTableWidgetItem(error);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);
    ui->sim_groupTbl->setItem(index, 4, item);

    ui->sim_groupTbl->selectRow(index);

}

