#include "simplotwindow.h"
#include "ui_simplotwindow.h"

SimPlotWindow::SimPlotWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SimPlotWindow)
{
    ui->setupUi(this);

    setWindowIcon(QIcon(":icon.svg"));
    connect(ui->simplot_plot1, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMovePlot1(QMouseEvent*)));
    connect(ui->simplot_plot1, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseClickPlot1(QMouseEvent*)));
    connect(ui->simplot_plot2, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMovePlot2(QMouseEvent*)));
    connect(ui->simplot_plot2, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseClickPlot2(QMouseEvent*)));

}

SimPlotWindow::~SimPlotWindow()
{
    delete ui;
}

void SimPlotWindow::initialize(std::vector<std::vector<Variant>> variants, std::vector<SimulationRequest> reqs, QString title){

    if(reqs[0].isRare())
        this->variants = filterCollapsed(variants, reqs[0].collapse);
    else
        this->variants = variants;
    this->requests = reqs;
    this->alpha = 0.05;
    this->ntests = variants[0][0].nPvals();
    for(int i = 0; i < ntests; i++)
        testTypes.push_back(QString::fromStdString(variants[0][0].getPvalSource(i)));

    if(reqs[0].underNull())
        this->yAxisLabel = "Type I Error";
    else
        this->yAxisLabel = "Power";

    this->xAxisLabel = "Test";
    this->setWindowTitle("Plotter - " + title);
    this->indexForPlot2 = -1;
    this->graphIndexForPlot2 = -1;

    QColor blue = QColor(106, 176, 203);
    QColor orange = QColor(255, 187, 63);
    QColor pink = QColor(204, 121, 167);
    QColor green = QColor(106, 150, 71);

    colours.push_back(blue);
    colours.push_back(orange);
    colours.push_back(pink);
    colours.push_back(green);

    buildPlot();
    buildLegend();
}

void SimPlotWindow::updateSampleSize(int index){

    if(index < 0){
        if(indexForPlot2 < 0){
            ui->simplot_ncontrolsDgt->display(0);
            ui->simplot_ncasesDgt->display(0);
            return;
        }
        else
            index = indexForPlot2;
    }

    ui->simplot_ncontrolsDgt->display(this->requests[index].ncase());
    ui->simplot_ncasesDgt->display(this->requests[index].ncontrol());

}

std::vector<std::vector<Variant>> SimPlotWindow::filterCollapsed(std::vector<std::vector<Variant>> variants, int k){
    std::vector<std::vector<Variant>> keep;

    for(int i =0; i < variants.size(); i++){
        std::vector<Variant> keep_i;

        int counter = 0;
        for(int j = 0; j < variants[i].size(); j++){

            counter++;

            if(counter == 1)
                keep_i.push_back(variants[i][j]);

            if(counter == k)
                counter = 0;
        }

        keep.push_back(keep_i);
    }

    return keep;
}

void SimPlotWindow::mouseClickPlot1(QMouseEvent *event){
    int closestIndex = findClosestPoint(ui->simplot_plot1, event);
    if(closestIndex >= 0){
        indexForPlot2 = closestIndex;
        buildPlot2(closestIndex);
    }
}

void SimPlotWindow::mouseClickPlot2(QMouseEvent *event){
    if(indexForPlot2 >= 0){
        int closestIndex = findClosestPoint(ui->simplot_plot2, event, true);
        graphIndexForPlot2 = closestIndex;
        buildPlot2(indexForPlot2, closestIndex);
    }
}

void SimPlotWindow::mouseMovePlot1(QMouseEvent *event){

    int closestIndex = findClosestPoint(ui->simplot_plot1, event);
    updateSampleSize(closestIndex);

    if(closestIndex >= 0){

        //remove old highlight layer
        ui->simplot_plot1->removeGraph(ui->simplot_plot1->graphCount()-1);
        ui->simplot_plot1->addGraph();

        QVector<double> x;
        QVector<double> y;

        for(int j = 0; j < ui->simplot_plot1->graphCount()-1; j++){
            x.push_back(ui->simplot_plot1->graph(j)->data().data()->at(closestIndex)->key);
            y.push_back(ui->simplot_plot1->graph(j)->data().data()->at(closestIndex)->value);
        }

        ui->simplot_plot1->graph()->setData(x, y);
        ui->simplot_plot1->graph()->setScatterStyle(
                    QCPScatterStyle(QCPScatterStyle::ssDisc,
                                    highlight, Qt::white, 8));
        ui->simplot_plot1->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);

        ui->simplot_plot1->replot();
    }
    else{
        //remove old highlight layer
        ui->simplot_plot1->removeGraph(ui->simplot_plot1->graphCount()-1);
        ui->simplot_plot1->addGraph();
        ui->simplot_plot1->replot();
    }
}

void SimPlotWindow::mouseMovePlot2(QMouseEvent *event){

    int closestGraphIndex = findClosestPoint(ui->simplot_plot2, event, true);
    buildPlot2(indexForPlot2, closestGraphIndex);
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


double SimPlotWindow::calculatePower(int index, std::vector<Variant> run, double alpha){
    double success = 0;

    for(Variant v : run)
        if(v.getPval(index) <= alpha)
            success++;

    return success/run.size();
}

QVector<double> SimPlotWindow::calculatePower(int index, std::vector<std::vector<Variant>> variants, double alpha){
    QVector<double> power;
    for(int i = 0; i < variants.size(); i++)
        power.push_back(calculatePower(index, variants[i], alpha));

    return power;
}

void SimPlotWindow::buildPlot(){

    ui->simplot_plot1->clearPlottables();
    ui->simplot_plot1->addGraph();

    QVector<double> index;
    for(int i = 1; i <= variants.size(); i++)
        index.push_back(i);

    for(int i = 0; i < ntests; i++){
        ui->simplot_plot1->graph()->setData(index, calculatePower(i, variants, alpha));
        ui->simplot_plot1->graph()->setPen(QPen(colours[i], 2));
        ui->simplot_plot1->graph()->setScatterStyle(QCPScatterStyle::ssDisc);
        //last layer is highlight layer
        ui->simplot_plot1->addGraph();
    }
    ui->simplot_plot1->rescaleAxes();
    ui->simplot_plot1->xAxis->setLabel(xAxisLabel);
    ui->simplot_plot1->yAxis->setLabel(yAxisLabel);
    ui->simplot_plot1->xAxis->setLabelFont(axisFont);
    ui->simplot_plot1->yAxis->setLabelFont(axisFont);
    ui->simplot_plot1->yAxis->setRange(0, 1.05);
    ui->simplot_plot1->xAxis->setRange(0.95, ui->simplot_plot1->xAxis->range().upper + 0.05);

    ui->simplot_plot1->replot();
}

void SimPlotWindow::buildLegend(){

    ui->simplot_legend->xAxis->setVisible(false);
    ui->simplot_legend->yAxis->setVisible(false);

    for(int i = 0; i < ntests; i++){
        ui->simplot_legend->addGraph();
        ui->simplot_legend->graph()->setPen(QPen(colours[i], 2));
        ui->simplot_legend->graph()->setScatterStyle(QCPScatterStyle::ssDisc);
        ui->simplot_legend->graph()->setName(testTypes[i]);
    }

    QCPLegend *legend = ui->simplot_legend->legend;
    legend->setVisible(true);
    legend->setFillOrder(QCPLayoutGrid::FillOrder::foColumnsFirst,true);
    legend->setFont(axisFont);
    legend->setBorderPen(Qt::NoPen);
    ui->simplot_legend->plotLayout()->addElement(1, 0, legend);
    ui->simplot_legend->plotLayout()->removeAt(0);
    ui->simplot_legend->repaint();

}


void SimPlotWindow::buildPlot2(int index, int focusGraph){

    if(index < 0)
        return;

    if(focusGraph < 0)
        focusGraph = graphIndexForPlot2;

    ui->simplot_plot2->clearItems();
    ui->simplot_plot2->clearPlottables();

    std::vector<Variant> toPlot = variants[index];

    double max = 0;

    for(int j = 0; j < ntests; j++){

        QVector<double> qq;
        QVector<double> p;

        for(int i = 0; i< toPlot.size(); i++){
            qq.push_back(-std::log10((i+1.0) * 1.0/(toPlot.size()*1.0)));
            p.push_back(-std::log10(toPlot[i].getPval(j)));
        }

        std::sort(p.begin(), p.end());
        std::reverse(qq.begin(), qq.end());

        ui->simplot_plot2->addGraph();
        ui->simplot_plot2->graph()->setData(qq, p);
        QPen scatterPen = QPen(colours[j], 2);

        if(focusGraph >= 0 && focusGraph != j)
            scatterPen = QPen(QColor(colours[j].red(), colours[j].green(), colours[j].blue(), 40), 2);

        ui->simplot_plot2->graph()->setPen(scatterPen);
        ui->simplot_plot2->graph()->setScatterStyle(QCPScatterStyle::ssDisc);
        ui->simplot_plot2->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);

        double testMax = std::max(p.last(), qq.last());
        if(testMax > max)
            max = testMax;
    }

    if(focusGraph >= 0){
        if(ui->simplot_plot2->layer("top") == ui->simplot_plot2->layer("dne"))
            ui->simplot_plot2->addLayer("top", ui->simplot_plot2->layer("main"),
                                        QCustomPlot::LayerInsertMode::limAbove);
        ui->simplot_plot2->graph(focusGraph)->setLayer("top");
    }

    //empty
    ui->simplot_plot2->addGraph();

    ui->simplot_plot2->xAxis->setLabel("-log10(Theoretical P-values)");
    ui->simplot_plot2->yAxis->setLabel("-log10(Observed P-values)");
    ui->simplot_plot2->xAxis->setLabelFont(axisFont);
    ui->simplot_plot2->yAxis->setLabelFont(axisFont);
    ui->simplot_plot2->rescaleAxes();

    ui->simplot_plot2->yAxis->setRange(0, max + max*0.01);
    ui->simplot_plot2->xAxis->setRange(0, max + max*0.01);

    QPen dpen = QPen(Qt::DashDotLine);
    dpen.setColor( QColor(0, 0, 0, 125));
    auto diagonal = new QCPItemLine(ui->simplot_plot2);
    diagonal->setPen(dpen);
    diagonal->start->setCoords(-1,-1);
    diagonal->end->setCoords(1e4,1e4);

    QPen hpen = QPen(Qt::DashDotLine);
    hpen.setColor(redLine);
    auto horizontal = new QCPItemLine(ui->simplot_plot2);
    horizontal->setPen(hpen);
    horizontal->start->setCoords(-1, -std::log10(alpha));
    horizontal->end->setCoords(1e4,-std::log10(alpha));

    ui->simplot_plot2->replot();
}

void SimPlotWindow::on_simplot_alphaDial_valueChanged(int value){
    alpha = std::pow(10, -value/1000.0);
    ui->simplot_alphaTxt->setText(QString::number(alpha));
    buildPlot();
    buildPlot2(indexForPlot2);
}

void SimPlotWindow::on_simplot_alphaTxt_textChanged(const QString &arg1){

    bool ok = false;
    double x = arg1.toDouble(&ok);
    if(ok && x >= 0 && x <=1){
        alpha = x;
        buildPlot();
        buildPlot2(indexForPlot2);
    }
}
