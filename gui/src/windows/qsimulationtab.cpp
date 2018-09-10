#include "MainWindow.h"
//#include "simplotwindow.h"

#include "ui_mainwindow.h"

const QString sep = ":";

std::vector<SimulationRequestGroup> MainWindow::qConstructGroups(int test, int ntest){


    std::vector<SimulationRequestGroup> groups;
    /*
    for(int i = 0; i< ui->sim_groupTbl->rowCount(); i++){

        SimulationRequestGroup g;

        QString size = ui->sim_groupTbl->item(i, 0)->text();
        int n;

        if(size.contains(sep)){

            bool ok1 = false;
            bool ok2 = false;

            int nMin = size.split(sep).at(0).toInt(&ok1);
            int nMax = size.split(sep).at(1).toInt(&ok2);
            if(!ok1 || !ok2 || nMin <= 0 || nMax <= 0)
                throwError("Invalid group settings. Ensure number of individuals is an integer greater than 0. If you wish to use a range please separate two integers by a colon \":\"");


            if(nMax < nMin){
                int temp = nMax;
                nMax = nMin;
                nMin = temp;
            }

            double denom = std::max(1, ntest-1);
            n = (nMax - nMin)/denom * test + nMin;
        }
        else
            n = size.toInt();

        QString caseControl = ((QComboBox*)ui->sim_groupTbl->indexWidget(ui->sim_groupTbl->model()->index(i, 1)))->currentText();
        QString meanDepth = ui->sim_groupTbl->item(i, 2)->text();
        QString sdDepth = ui->sim_groupTbl->item(i, 3)->text();
        QString errorRate = ui->sim_groupTbl->item(i, 4)->text();

        g.n = n;
        g.index = i;

        bool ok = false;

        g.meanDepth = meanDepth.toDouble(&ok);
        if( !ok || g.meanDepth <= 0)
            throwError("Invalid group settings. Expected mean read depth to be a numeric value greater than 0.");


        ok=false;
        g.sdDepth = sdDepth.toDouble(&ok);
        if(!ok || g.sdDepth < 0)
            throwError("Invalid group settings. Expected read depth SD to be a numeric value greater than or equal to 0.");

        ok=false;
        g.errorRate = errorRate.toDouble(&ok);
        if(!ok || g.errorRate < 0 || g.errorRate > 1)
            throwError("Invalid group settings. Expected error rate to be a numeric value between 0 and 1");


        g.isCase = (caseControl == "case");

        ok=false;
        int highLow = ui->sim_groupHighLowTxt->text().toInt(&ok);

        if(!ok || highLow <= 0)
            throwError("Expected high-low read depth cutoff to be an integer greater than 0");


        g.isHrg = g.meanDepth < highLow;
        groups.push_back(g);
    }

    */
    return groups;
}

SimulationRequest MainWindow::qConstructRequest(std::vector<SimulationRequestGroup> groups){

    SimulationRequest request;

    /*
    request.groups = groups;

    request.test = "common";
    request.useBootstrap = false;
\
    if(ui->sim_testRareCalphaBtn->isChecked()){
        request.test = "calpha";
        request.collapse = ui->sim_collapseTxt->text().toInt();
        request.useBootstrap = true;

    }
    else if(ui->sim_testRareCastBtn->isChecked()){
        request.test= "cast";
        request.collapse = ui->sim_collapseTxt->text().toInt();
        request.useBootstrap = true;
    }

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

    */
    return request;
}


void MainWindow::on_qsim_stopBtn_clicked(){
    stopJob();
}

void MainWindow::on_qsim_runBtn_clicked(){
    if(!ui->qsim_runBtn->isEnabled())
        return;

    /*
    greyOutput();
    disableRun();

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
        jobThread = thread;
    }catch(...){
        ui->sim_runBtn->setEnabled(true);
        enableRun();
    }
    */
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
    ui->sim_collapseTxt->setEnabled(value);
    ui->sim_collapseLbl->setEnabled(value);
    ui->sim_testBootLbl->setEnabled(value);
    ui->sim_testBootTxt->setEnabled(value);
    ui->sim_testStopChk->setEnabled(value);
}

void MainWindow::on_qsim_testRareCastBtn_toggled(bool checked){
    qSimEnableRare(checked);
}

void MainWindow::on_qsim_testRareSkatBtn_toggled(bool checked){
    qSimEnableRare(checked);
}

void MainWindow::qAddGroup(QString n, bool control, QString depth, QString sdDepth, QString errorRate){

    /*
    ui->sim_groupTbl->setSelectionBehavior(QAbstractItemView::SelectRows);

    int index = ui->sim_groupTbl->rowCount();

    ui->sim_groupTbl->insertRow(index);

    QTableWidgetItem* item;

    item = new QTableWidgetItem(n);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsEditable);
    ui->sim_groupTbl->setItem(index, 0, item);

    QComboBox *caseControlBox = new QComboBox();
    caseControlBox->addItem("case");
    caseControlBox->addItem("control");

    if(control)
        caseControlBox->setCurrentIndex(1);

    ui->sim_groupTbl->setIndexWidget(ui->sim_groupTbl->model()->index(index, 1), caseControlBox);

    item = new QTableWidgetItem(depth);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsEditable);
    ui->sim_groupTbl->setItem(index, 2, item);

    item = new QTableWidgetItem(sdDepth);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsEditable);
    ui->sim_groupTbl->setItem(index, 3, item);

    item = new QTableWidgetItem(errorRate);
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsEditable);
    ui->sim_groupTbl->setItem(index, 4, item);

    ui->sim_groupTbl->selectRow(index);
    ui->sim_groupTbl->update();
    */
}

void MainWindow::on_qsim_groupAddBtn_clicked(){
    //addGroup("100", true, "10", "2", "0.01");
}

