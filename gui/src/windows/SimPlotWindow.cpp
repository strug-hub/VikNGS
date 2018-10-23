#include "SimPlotWindow.h"
#include "ui_simplotwindow.h"
#include "TableDisplayWindow.h"

SimPlotWindow::SimPlotWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SimPlotWindow)
{
    ui->setupUi(this);

    setWindowIcon(QIcon(":icon.svg"));
    connect(ui->simplot_power, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMovePlot1(QMouseEvent*)));
    connect(ui->simplot_power, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseClickPlot1(QMouseEvent*)));
    connect(ui->simplot_plot2, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMovePlot2(QMouseEvent*)));
    connect(ui->simplot_plot2, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseClickPlot2(QMouseEvent*)));

    QColor blue = QColor(106, 176, 203);
    QColor orange = QColor(255, 187, 63);
    QColor pink = QColor(204, 121, 167);
    QColor green = QColor(106, 150, 71);

    colours.push_back(blue);
    colours.push_back(pink);
    colours.push_back(orange);
    colours.push_back(green);
}

SimPlotWindow::~SimPlotWindow(){
    delete ui;
}

void SimPlotWindow::initialize(Data& results, SimulationRequest& req, QString title){

    this->nsteps = req.steps;
    this->testsPerStep = std::round(1.0 * results.tests.size() / (1.0 * nsteps));
    this->alpha = 0.05;

    if(req.underNull()){
        this->yAxisLabel = "Type I Error";
        ui->simplot_measureGrp->setTitle("Type I Error");
    }
    else{
        this->yAxisLabel = "Power";
        ui->simplot_measureGrp->setTitle("Power");
    }

    this->xAxisLabel = "Test";
    this->setWindowTitle("Plotter - " + title);

    getPvalues(results.variants);

    this->result = results;
    this->request = req;

    this->powerIndex = 0;
    this->qqIndex = -1;

    buildPowerPlot();
    buildLegend();

    buildQQPlot(powerIndex);
    updatePowerValues(0);
}

void SimPlotWindow::getPvalues(std::vector<VariantSet>& variants){
    if(variants.size() < 1)
        return;

    for(size_t i = 0; i < variants[0].nPvals(); i++){
        std::vector<double> pval_i(variants.size());

        for(size_t j = 0; j < variants.size(); j++){
            if(variants[j].nPvals() > i)
                pval_i[j] = variants[j].getPval(i);
        }
        std::sort(pval_i.begin(), pval_i.end());
        pvalues.push_back(pval_i);
    }
}

void SimPlotWindow::updateSampleSize(int index){

    if(index < 0){
        if(powerIndex < 0){
            ui->simplot_ncontrolsDgt->display(0);
            ui->simplot_ncasesDgt->display(0);
            return;
        }
        else
            index = powerIndex;
    }

    ui->simplot_ncasesDgt->display(this->request.ncase(powerIndex));
    ui->simplot_ncontrolsDgt->display(this->request.ncontrol(powerIndex));

}

void SimPlotWindow::mouseClickPlot1(QMouseEvent *event){
    int closestIndex = findClosestPoint(ui->simplot_power, event);
    if(closestIndex >= 0){
        powerIndex = closestIndex;
        buildQQPlot(closestIndex);
        updatePowerValues(closestIndex);
    }
}

void SimPlotWindow::mouseClickPlot2(QMouseEvent *event){
    if(powerIndex >= 0){
        int closestIndex = findClosestPoint(ui->simplot_plot2, event, true);
        qqIndex = closestIndex;
        buildQQPlot(powerIndex, closestIndex);
    }
}

void SimPlotWindow::mouseMovePlot1(QMouseEvent *event){

    int closestIndex = findClosestPoint(ui->simplot_power, event);
    updateSampleSize(closestIndex);

    if(closestIndex >= 0){

        //remove old highlight layer
        ui->simplot_power->removeGraph(ui->simplot_power->graphCount()-1);
        ui->simplot_power->addGraph();

        QVector<double> x;
        QVector<double> y;

        for(int j = 0; j < ui->simplot_power->graphCount()-1; j++){
            x.push_back(ui->simplot_power->graph(j)->data().data()->at(closestIndex)->key);
            y.push_back(ui->simplot_power->graph(j)->data().data()->at(closestIndex)->value);
        }

        ui->simplot_power->graph()->setData(x, y);
        ui->simplot_power->graph()->setScatterStyle(
                    QCPScatterStyle(QCPScatterStyle::ssDisc,
                                    highlight, Qt::white, 8));
        ui->simplot_power->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);

        ui->simplot_power->replot();
    }
    else{
        //remove old highlight layer
        ui->simplot_power->removeGraph(ui->simplot_power->graphCount()-1);
        ui->simplot_power->addGraph();
        ui->simplot_power->replot();
    }
}

void SimPlotWindow::mouseMovePlot2(QMouseEvent *event){

    int closestGraphIndex = findClosestPoint(ui->simplot_plot2, event, true);

    if(qqIndex != closestGraphIndex)
        buildQQPlot(powerIndex, closestGraphIndex);
}

int SimPlotWindow::findClosestPoint(QCustomPlot *plot, QMouseEvent *event, bool getGraphIndex){

    double maxDist = 200;

    double x = event->pos().x();
    double y = event->pos().y();

    double minDist = 9999999999;
    int minIndex;

    for(int j = 0; j < plot->graphCount()-1; j++){
        for(int i = 0; i <  plot->graph(j)->dataCount(); i++){

            double dataCoordx = plot->graph(j)->data().data()->at(i)->key;
            double dataCoordy = plot->graph(j)->data().data()->at(i)->value;

            double dataX = plot->xAxis->coordToPixel(dataCoordx);
            double dataY = plot->yAxis->coordToPixel(dataCoordy);

            double dist = pow(dataX - x, 2) + pow(dataY - y, 2);

            if(dist < minDist){
                minDist = dist;
                if(getGraphIndex)
                    minIndex = j;
                else
                    minIndex = i;
            }
        }
    }
    if(minDist > maxDist)
        return -1;

    return minIndex;
}

double SimPlotWindow::calculatePower(int index, double alpha){
    double success = 0;
    double npval = pvalues[index].size();
    for(size_t i = 0; i < npval; i++)
        if(pvalues[index][i] <= alpha)
            success++;

    return success/npval;
}

void SimPlotWindow::updatePowerValues(int index){

    ui->simplot_powerTbl->setColumnCount(2);
    ui->simplot_powerTbl->setRowCount(testsPerStep);

    for(int i = 0; i < testsPerStep; i++){
        double power = ui->simplot_power->graph(i)->data().data()->at(index)->value;
        ui->simplot_powerTbl->setItem(i, 1, new QTableWidgetItem(QString::number(power)));
        QTableWidgetItem * c = new QTableWidgetItem(QString::fromStdString(result.tests[i].toShortString()));
        c->setBackgroundColor(colours[i]);
        ui->simplot_powerTbl->setItem(i, 0, c);
    }
}

void SimPlotWindow::on_simplot_alphaDial_valueChanged(int value){
    alpha = std::pow(10, -value/1000.0);
    ui->simplot_alphaTxt->setText(QString::number(alpha));
    buildPowerPlot();
    updateAlphaLine();
    updatePowerValues(powerIndex);
}

void SimPlotWindow::on_simplot_alphaTxt_textChanged(const QString &arg1){

    bool ok = false;
    double x = arg1.toDouble(&ok);
    if(ok && x >= 0 && x <=1){
        alpha = x;

        buildPowerPlot();
        updateAlphaLine();
        updatePowerValues(powerIndex);
    }
}

void SimPlotWindow::on_pushButton_clicked(bool checked)
{
    TableDisplayWindow *table = new TableDisplayWindow();
    QString title = "Table";

    std::vector<int> testsToShow;

    for(int i = 0; i < this->testsPerStep; i++)
        testsToShow.push_back(this->powerIndex * this->testsPerStep + i);

    table->initialize(title, &result, testsToShow);
    table->show();
}
