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

            n = (nMax - nMin)/ntest * test + nMin;
        }
        else
            n = size.toInt();

        QString caseControl = ui->sim_groupTbl->item(i, 1)->text();
        QString meanDepth = ui->sim_groupTbl->item(i, 2)->text();
        QString sdDepth = ui->sim_groupTbl->item(i, 3)->text();

        QString highLow =
                ui->sim_groupTbl->item(i, 2)->text().toDouble() < ui->sim_groupHighLowTxt->text().toDouble() ?
                    "low" : "high";

        g.index = i;
        g.n = n;
        g.isCase = (caseControl == "case");
        g.isHrg = (highLow == "high");
        g.meanDepth = meanDepth.toDouble();
        g.sdDepth = sdDepth.toDouble();

        groups.push_back(g);

    }

    return groups;
}

SimulationRequest MainWindow::constructRequest(std::vector<SimulationRequestGroup> groups){

    SimulationRequest request;
    request.groups = groups;

    request.test = "common";
    if(ui->sim_testRareCalphaBtn->isChecked())
        request.test = "calpha";
    else if(ui->sim_testRareCastBtn->isChecked())
        request.test= "cast";

    request.useBootstrap = ui->sim_testBootChk->isChecked();
    request.stopEarly = false;
    if (request.useBootstrap){
        request.nboot = ui->sim_testBootTxt->text().toInt();
        request.stopEarly = ui->sim_testStopChk->isChecked();
    }
    else
        request.nboot=0;

    request.npop = ui->sim_populationSizeTxt->text().toInt();
    request.prevalence = ui->sim_populationPrevalenceTxt->text().toDouble();

    request.nsnp = ui->sim_variantSizeTxt->text().toInt();
    request.me = ui->sim_sequencingErrorRateTxt->text().toDouble();
    request.sde = ui->sim_sequencingErrorSdTxt->text().toDouble();

    request.oddsRatio = ui->sim_oddsRatioTxt->text().toDouble();
    request.maf = ui->sim_variantMafTxt->text().toDouble();

    return request;
}

double MainWindow::calculatePower( QVector<double> pval){

    int total = pval.size();
    int count = 0;
    double alpha = ui->sim_powerAlphaTxt->text().toDouble();

    for(int i = 0; i < total; i++)
        if(pval.at(i) <= alpha)
            count++;

    std::cout << (1.0*count)/(1.0*total);
    std::cout << "\n";

    return (1.0*count)/(1.0*total);
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
}

void MainWindow::on_sim_testRareCalphaBtn_toggled(bool checked){
    if(checked)
        ui->sim_testBootChk->setChecked(true);
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

    bool ok = false;
    double checkMean = mean.toDouble(&ok);

    if( !ok || checkMean <= 0){
        warningDialog("Invalid group settings. Expected mean read depth to be a numeric value greater than 0.");
        return;
    }

    ok=false;
    int checkSd = sd.toInt(&ok);

    if(!ok || checkSd < 0){
        warningDialog("Invalid group settings. Expected read depth SD to be a numeric value greater than or equal to 0.");
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


    ui->sim_groupTbl->selectRow(index);

}

