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

}

SimPlotWindow::~SimPlotWindow()
{
    delete ui;
}

void SimPlotWindow::mouseClickPlot1(QMouseEvent *event){
    int closestIndex = findClosestPoint(ui->simplot_plot1, event);
    if(closestIndex >= 0)
        buildPlot2(closestIndex);
}

void SimPlotWindow::buildPlot2(int index){

    std::vector<Variant> toPlot = variants[index];

    QVector<double> qq;
    QVector<double> p;

    for(int i = 0; i< toPlot.size(); i++){
        qq.push_back((i+1.0) * 1.0/(toPlot.size()*1.0));
        p.push_back(toPlot[i].pvalue);
    }

    std::sort(p.begin(), p.end());

    ui->simplot_plot2->clearPlottables();
    ui->simplot_plot2->addGraph();
    ui->simplot_plot2->xAxis->setLabel("Observed");
    ui->simplot_plot2->yAxis->setLabel("Theoretical");
    ui->simplot_plot2->graph()->setData(p, qq);
    ui->simplot_plot2->graph()->setScatterStyle(QCPScatterStyle::ssDisc);
    ui->simplot_plot2->graph()->setLineStyle( QCPGraph::LineStyle::lsNone);

    ui->simplot_plot2->rescaleAxes();
    ui->simplot_plot2->yAxis->setRange(0, 1);
    ui->simplot_plot2->xAxis->setRange(0, 1);

    //highlight layer
    ui->simplot_plot2->addGraph();

    QPen dpen = QPen(Qt::DashDotLine);
    dpen.setColor( QColor::fromRgb(210, 80, 80));
    auto diagonal = new QCPItemLine(ui->simplot_plot2);
    diagonal->setPen(dpen);
    diagonal->start->setCoords(-1,-1);
    diagonal->end->setCoords(100,100);

    ui->simplot_plot2->replot();
}

void SimPlotWindow::mouseMovePlot1(QMouseEvent *event){

    int closestIndex = findClosestPoint(ui->simplot_plot1, event);
    if(closestIndex >= 0){

        //remove old highlight layer
        ui->simplot_plot1->removeGraph(ui->simplot_plot1->graphCount()-1);
        ui->simplot_plot1->addGraph();

        QVector<double> x;
        QVector<double> y;

        x.push_back(ui->simplot_plot1->graph(0)->data().data()->at(closestIndex)->key);
        y.push_back(ui->simplot_plot1->graph(0)->data().data()->at(closestIndex)->value);

        ui->simplot_plot1->graph()->setData(x, y);
        ui->simplot_plot1->graph()->setScatterStyle(
                    QCPScatterStyle(QCPScatterStyle::ssDisc,
                                    highlight, Qt::white, 8));
        ui->simplot_plot1->replot();
    }
    else{
        //remove old highlight layer
        ui->simplot_plot1->removeGraph(ui->simplot_plot1->graphCount()-1);
        ui->simplot_plot1->addGraph();
        ui->simplot_plot1->replot();
    }

}

int SimPlotWindow::findClosestPoint(QCustomPlot *plot, QMouseEvent *event){

    double maxDist = 200;

    double x = event->pos().x();
    double y = event->pos().y();

    double minDist = 9999999999;
    int minIndex;

    int npoints = plot->graph(0)->data().data()->dataRange().end();
    for(int i = 0; i<npoints; i++){

        double dataCoordx = plot->graph(0)->data().data()->at(i)->key;
        double dataCoordy = plot->graph(0)->data().data()->at(i)->value;

        double dataX = plot->xAxis->coordToPixel(dataCoordx);
        double dataY = plot->yAxis->coordToPixel(dataCoordy);

        double dist = pow(dataX - x, 2) + pow(dataY - y, 2);

        if(dist < minDist){
            minDist = dist;
            minIndex = i;
        }
    }

    if(minDist > maxDist)
        return -1;

    return minIndex;
}


void SimPlotWindow::initialize(std::vector<std::vector<Variant>> variants, std::vector<SimulationRequest> reqs, QString title){

    this->variants = variants;
    this->requests = reqs;
    this->setWindowTitle("Plotter - " + title);
    buildPlot();
}

double SimPlotWindow::calculatePower(std::vector<Variant> run, double alpha){
    double success = 0;

    for(Variant v : run)
        if(v.pvalue <= alpha)
            success++;

    return success/run.size();
}

QVector<double> SimPlotWindow::calculatePower(std::vector<std::vector<Variant>> variants, double alpha){
    QVector<double> power;
    for(int i = 0; i< variants.size(); i++)
        power.push_back(calculatePower(variants[i], alpha));

    return power;
}

void SimPlotWindow::buildPlot(){
    double alpha = 0.05;

    ui->simplot_plot1->clearPlottables();
    ui->simplot_plot1->addGraph();
    ui->simplot_plot1->xAxis->setLabel("Test");
    ui->simplot_plot1->yAxis->setLabel("Power");

    QVector<double> index;
    for(int i = 1; i <= variants.size(); i++)
        index.push_back(i);

    ui->simplot_plot1->graph()->setData(index, calculatePower(variants, alpha));
    ui->simplot_plot1->graph()->setScatterStyle(QCPScatterStyle::ssDisc);

    ui->simplot_plot1->rescaleAxes();
    ui->simplot_plot1->yAxis->setRange(0, 1);

    //highlight layer
    ui->simplot_plot1->addGraph();

    ui->simplot_plot1->replot();


    /*
    focusedChr = chrName;
    Chromosome chr = chromosomes[chrName];

    ui->plot_chrPlt->addGraph();

    QVector<double> pos = chr.getPositions();
    QVector<double> pvals = chr.getPvals();

    QVector<double> pos_range;
    QVector<double> pvals_range;

    for(int i =0; i<pos.size(); i++){

        if(min < 0 || pos[i] > min)
            if(max < 0 || pos[i] < max){
                pos_range.push_back(pos[i]);
                pvals_range.push_back(pvals[i]);
            }
    }

    ui->plot_chrPlt->graph()->setData(pos_range,pvals_range);
    ui->plot_chrPlt->graph()->setName(chrGraphName);
    ui->plot_chrPlt->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);

    double pointSize = 1.0/(sqrt(pos_range.size())/100.0);
    if(pointSize > 6) pointSize = 6;
    if(pointSize < 2) pointSize = 2;

    ui->plot_chrPlt->graph()->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                focus, Qt::white, pointSize));

    ui->plot_chrPlt->yAxis->setLabel("-log(p)");
    ui->plot_chrPlt->xAxis->setLabel("Position");

    QCPTextElement *title = new QCPTextElement(ui->plot_chrPlt, "Chromosome " + chrName, QFont("sans", 17, QFont::Bold));
    ui->plot_chrPlt->plotLayout()->removeAt(0);
    ui->plot_chrPlt->plotLayout()->addElement(0,0,title);

    ui->plot_chrPlt->rescaleAxes();
    if(ui->plot_chrPlt->yAxis->atTop < 7.5)
        ui->plot_chrPlt->yAxis->setRange(0, 7.5);

    QPen hpen = QPen(Qt::DashDotLine);
    hpen.setColor( QColor::fromRgb(210, 80, 80));
    horizontal = new QCPItemLine(ui->plot_chrPlt);
    horizontal->setPen(hpen);
    horizontal->start->setCoords(0,1.301);
    horizontal->end->setCoords(chr.getMaxPos(),1.301);

    QPen bpen = QPen(Qt::DashDotLine);
    bpen.setColor( QColor::fromRgb(250, 145, 145));
    QCPItemLine *bonferroni = new QCPItemLine(ui->plot_chrPlt);
    bonferroni->setPen(bpen);
    bonferroni->start->setCoords(0,7.301);
    bonferroni->end->setCoords(chr.getMaxPos(), 7.301);

    ui->plot_chrPlt->addGraph();
    ui->plot_chrPlt->graph()->setName("annotationLayer");
    ui->plot_chrPlt->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);
    ui->plot_chrPlt->graph()->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                highlight, Qt::white, 6));

    updateVariantHighlightLayer(nullVariant);

    if(min < 0)
        ui->plot_zoomBar->updateMin(ui->plot_chrPlt->xAxis->range().lower);
    else
        ui->plot_zoomBar->updateZoomMin(min);

    if(max < 0)
        ui->plot_zoomBar->updateMax(ui->plot_chrPlt->xAxis->range().upper);
    else
        ui->plot_zoomBar->updateZoomMax(max);


     ui->plot_chrPlt->replot();
     */
}



