#include "mainwindow.h"
#include "ui_mainwindow.h"

std::vector<SimulationRequestGroup> MainWindow::constructGroups(int test, int ntest){

    std::vector<SimulationRequestGroup> groups;

    try{
        for(int i = 0; i< ui->sim_groupTbl->rowCount(); i++){

            int groupID = i;
            QString size = ui->sim_groupTbl->item(i, 0)->text();
            std::string n;

            if(size.contains(':')){
                int nMin = size.split(':').at(0).toInt();
                int nMax = size.split(':').at(1).toInt();

                n = std::to_string((nMax - nMin)/ntest * test + nMin);
            }
            else
                n = size.toStdString();

            std::string caseControl = ui->sim_groupTbl->item(i, 1)->text().toStdString();
            std::string meanDepth = ui->sim_groupTbl->item(i, 2)->text().toStdString();
            std::string sdDepth = ui->sim_groupTbl->item(i, 3)->text().toStdString();

            std::string highLow =
                    ui->sim_groupTbl->item(i, 2)->text().toDouble() < ui->sim_groupHighLowTxt->text().toDouble() ?
                        "low" : "high";


            groups.push_back(
                        newSimulationRequestGroup( groupID,
                            n, caseControl, highLow, meanDepth, sdDepth)
                        );

            groups[i].print();
        }
    }
    catch(...){
        throw;
    }

    return groups;
}

SimulationRequest MainWindow::constructRequest(std::vector<SimulationRequestGroup> groups){

    SimulationRequest request;

    try{
        std::string test = "common";
        if(ui->sim_testRareCalphaBtn->isChecked())
            test = "calpha";
        else if(ui->sim_testRareCastBtn->isChecked())
            test= "cast";

        request = newSimulationRequest(
                    ui->sim_populationSizeTxt->text().toStdString(),
                    ui->sim_populationPrevalenceTxt->text().toStdString(),
                    ui->sim_variantSizeTxt->text().toStdString(),
                    ui->sim_sequencingErrorRateTxt->text().toStdString(),
                    ui->sim_sequencingErrorSdTxt->text().toStdString(),
                    ui->sim_oddsRatioTxt->text().toStdString(),
                    ui->sim_variantMafTxt->text().toStdString(),
                    groups,
                    test,
                    ui->sim_testBootChk->isChecked(),
                    ui->sim_testBootTxt->text().toStdString()
                    );

        request.print();
    }
    catch(...){
        throw;
    }

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


void MainWindow::on_sim_runBtn_clicked(){

    if(!ui->sim_runBtn->isEnabled())
        return;

    ui->sim_runBtn->setEnabled(false);

    int ntest = 1;
    bool multitest = false;
    for(int i = 0; i< ui->sim_groupTbl->rowCount(); i++)
        if(ui->sim_groupTbl->item(i, 0)->text().contains(':'))
            multitest =true;

    if(multitest)
        ntest = std::max(1, ui->sim_powerStepTxt->text().toInt());


    std::cout<< "RUN BABY RUN\n\n------\n";

    QVector<double> power;

    try{

        for(int run = 0; run < ntest; run++){

            std::vector<SimulationRequestGroup> groups = constructGroups(run, ntest);
            SimulationRequest request = constructRequest(groups);

            QFuture<std::vector<double>> future = QtConcurrent::run(startSimulation, request);

            while(future.isRunning()){  }

            QVector<double> pval = QVector<double>::fromStdVector(future.result());

            power.push_back(calculatePower(pval));
            sim_replot(power, ntest);
        }

    }catch(std::invalid_argument e){
        QMessageBox *dialog = new QMessageBox;
        dialog->setWindowTitle("Input error");

        QString error(e.what());

        dialog->setText("Invalid input. Expected " +
                        error.split(",").at(1) +
                        " for following value: " +
                        error.split(",").at(0) + ".");
        dialog->show();

        ui->sim_runBtn->setEnabled(true);

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

    ui->sim_runBtn->setEnabled(true);

}

void MainWindow::on_sim_groupAddBtn_clicked()
{
    ui->sim_groupTbl->setSelectionBehavior(QAbstractItemView::SelectRows);

    bool isCase = ui->sim_groupCaseBtn->isChecked();

    QString nMin = ui->sim_groupSizeMinTxt->text();
    QString nMax = ui->sim_groupSizeMaxTxt->text();
    QString n = (nMin == nMax) ? nMin : nMin + ":" + nMax;

    QString mean = ui->sim_groupDepthMeanTxt->text();
    QString sd = ui->sim_groupDepthSdTxt->text();

    int mn = nMin.toInt();
    int mx = nMax.toInt();

    if(mn <= 0 || mx <= 0){
        QMessageBox *dialog = new QMessageBox;
        dialog->setText("Invalid group settings. Ensure min and max size are integers greater than 0.");
        dialog->show();
        return;
    }

    if(mn > mx){
        QMessageBox *dialog = new QMessageBox;
        dialog->setText("Invalid group settings. Expected min size to be less than or equal to max size.");
        dialog->show();
        return;
    }

    if( mean.toDouble() <= 0){
        QMessageBox *dialog = new QMessageBox;
        dialog->setText("Invalid group settings. Expected mean read depth to be a numeric value greater than 0");
        dialog->show();
        return;
    }

   bool ok;
   int checkSd = sd.toDouble(&ok);

   if(!ok || checkSd < 0){
       QMessageBox *dialog = new QMessageBox;
       dialog->setText("Invalid group settings. Expected read depth SD to be a numeric value greater than or equal to 0");
       dialog->show();
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

void MainWindow::sim_replot(QVector<double> values, int maxX){

    int n = values.size();

    if(n < 0)
        return;

    QVector<double> x(n), y(n);

    for (int i=0; i < n; ++i)
    {
      x[i] = i;
      y[i] = values.at(i);
    }

    ui->sim_plotbox->clearPlottables();
    ui->sim_plotbox->addGraph();

    // give the axes some labels:
    ui->sim_plotbox->xAxis->setLabel("Test");
    ui->sim_plotbox->yAxis->setLabel("Power");

    ui->sim_plotbox->graph()->setData(x, y);
    ui->sim_plotbox->graph()->setScatterStyle(QCPScatterStyle::ssDisc);
    // set axes ranges, so we see all data:
    ui->sim_plotbox->xAxis->setRange(0, maxX);
    ui->sim_plotbox->yAxis->setRange(0, 1);
    ui->sim_plotbox->replot();
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

