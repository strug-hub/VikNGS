#include "simplotwindow.h"
#include "ui_simplotwindow.h"
#include "tabledisplaywindow.h"

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

void SimPlotWindow::getPvalues(std::vector<std::vector<Variant>> &variants){

    for(int k = 0; k < ntests; k++){

        std::vector<std::vector<double>> pval_i;

        for(int i = 0; i < nsteps; i++){
            std::vector<double> pval_ij;

            for(int j = 0; j < variants[i].size(); j++)
                pval_ij.push_back(variants[i][j].getPval(k));

            std::sort(pval_ij.begin(), pval_ij.end());
            pval_i.push_back(pval_ij);
        }
        pvalues.push_back(pval_i);
    }
    return;
}

void SimPlotWindow::initialize(std::vector<std::vector<Variant>>& variants, SimulationRequest& req, QString title){

    this->ntests = variants[0][0].nPvals();
    this->nsteps = variants.size();

    if(req.isRare())
        variants = filterCollapsed(variants, req.collapse);

    getPvalues(variants);
    this->variants = variants;

    this->request = req;
    this->alpha = 0.05;
    for(int i = 0; i < ntests; i++)
        testTypes.push_back(QString::fromStdString(variants[0][0].getPvalSource(i)));

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
    this->stepIndexForPlot2 = -1;
    this->testIndexForPlot2 = -1;

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

    stepIndexForPlot2 = 0;
    buildPlot2(stepIndexForPlot2);
    updatePowerValues(0);

}

void SimPlotWindow::updateSampleSize(int index){

    if(index < 0){
        if(stepIndexForPlot2 < 0){
            ui->simplot_ncontrolsDgt->display(0);
            ui->simplot_ncasesDgt->display(0);
            return;
        }
        else
            index = stepIndexForPlot2;
    }

    ui->simplot_ncasesDgt->display(this->request.ncase(stepIndexForPlot2));
    ui->simplot_ncontrolsDgt->display(this->request.ncontrol(stepIndexForPlot2));

}

std::vector<std::vector<Variant>> SimPlotWindow::filterCollapsed(std::vector<std::vector<Variant>> &variants, int k){
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

void SimPlotWindow::on_pushButton_pressed()
{

    TableDisplayWindow *table = new TableDisplayWindow();
    QString title = "Table";
    table->initialize(title, variants[stepIndexForPlot2], request, stepIndexForPlot2);
    table->show();
}

void SimPlotWindow::mouseClickPlot1(QMouseEvent *event){
    int closestIndex = findClosestPoint(ui->simplot_plot1, event);
    if(closestIndex >= 0){
        stepIndexForPlot2 = closestIndex;
        buildPlot2(closestIndex);
        updatePowerValues(closestIndex);
    }
}

void SimPlotWindow::mouseClickPlot2(QMouseEvent *event){
    if(stepIndexForPlot2 >= 0){
        int closestIndex = findClosestPoint(ui->simplot_plot2, event, true);
        testIndexForPlot2 = closestIndex;
        buildPlot2(stepIndexForPlot2, closestIndex);
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

    if(testIndexForPlot2 != closestGraphIndex)
        buildPlot2(stepIndexForPlot2, closestGraphIndex);
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

QVector<double> SimPlotWindow::calculatePower(int testIndex, double alpha){
    QVector<double> power;
    for(int i = 0; i < nsteps; i++){

        double success = 0;
        for(int j = 0; j < pvalues[testIndex][i].size(); j++)
            if(pvalues[testIndex][i][j] <= alpha)
                success++;
            else
                break;

        power.push_back(success/pvalues[testIndex][i].size());
    }

    return power;
}

void SimPlotWindow::updatePowerValues(int index){

    ui->simplot_powerTbl->setColumnCount(2);
    ui->simplot_powerTbl->setRowCount(ntests);

    for(int i = 0; i < ntests; i++){
        double power = ui->simplot_plot1->graph(i)->data().data()->at(index)->value;
        ui->simplot_powerTbl->setItem(i, 1, new QTableWidgetItem(QString::number(power)));
        QTableWidgetItem * c = new QTableWidgetItem(QString::fromStdString(variants[index][0].getPvalSourceShort(i)));
        c->setBackgroundColor(colours[i]);
        ui->simplot_powerTbl->setItem(i, 0, c);
    }
}


void SimPlotWindow::buildPlot(){

    ui->simplot_plot1->clearPlottables();
    ui->simplot_plot1->addGraph();

    QVector<double> index;
    for(int i = 1; i <= nsteps; i++)
            index.push_back(i);

    for(int i = 0; i < ntests; i++){
            QVector<double> power = calculatePower(i, alpha);
            ui->simplot_plot1->graph()->setData(index, power);
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

    if(ntests > 0){
        ui->simplot_legend1Lbl->setText(testTypes[0]);
        QPalette pal = palette();
        pal.setColor(QPalette::Background, colours[0]);
        ui->simplot_legend1Sqr->setAutoFillBackground(true);
        ui->simplot_legend1Sqr->setPalette(pal);
    }
    else{
        ui->simplot_legend1Lbl->hide();
        ui->simplot_legend1Sqr->hide();
    }

    if(ntests > 1){
        ui->simplot_legend2Lbl->setText(testTypes[1]);
        QPalette pal = palette();
        pal.setColor(QPalette::Background, colours[1]);
        ui->simplot_legend2Sqr->setAutoFillBackground(true);
        ui->simplot_legend2Sqr->setPalette(pal);
    }
    else{
        ui->simplot_legend2Lbl->hide();
        ui->simplot_legend2Sqr->hide();
    }

    if(ntests > 2){
        ui->simplot_legend3Lbl->setText(testTypes[2]);
        QPalette pal = palette();
        pal.setColor(QPalette::Background, colours[2]);
        ui->simplot_legend3Sqr->setAutoFillBackground(true);
        ui->simplot_legend3Sqr->setPalette(pal);
    }
    else{
        ui->simplot_legend3Lbl->hide();
        ui->simplot_legend3Sqr->hide();
    }

    if(ntests > 3){
        ui->simplot_legend4Lbl->setText(testTypes[3]);
        QPalette pal = palette();
        pal.setColor(QPalette::Background, colours[3]);
        ui->simplot_legend4Sqr->setAutoFillBackground(true);
        ui->simplot_legend4Sqr->setPalette(pal);
    }
    else{
        ui->simplot_legend4Lbl->hide();
        ui->simplot_legend4Sqr->hide();
    }
}

void SimPlotWindow::updateAlphaLine(){

    if(stepIndexForPlot2 >= 0){
        alphaLine->start->setCoords(-1, -std::log10(alpha));
        alphaLine->end->setCoords(1e4,-std::log10(alpha));

        ui->simplot_plot2->layer("overlay")->replot();
    }
}

void SimPlotWindow::buildPlot2(int stepIndex, int focusGraph){

    if(stepIndex < 0)
        return;

    if(focusGraph < 0)
        focusGraph = testIndexForPlot2;

    ui->simplot_plot2->clearItems();
    ui->simplot_plot2->clearPlottables();

    double max = 0;

    for(int testIndex = 0; testIndex < ntests; testIndex++){

        QVector<double> qq;
        QVector<double> p;

        double step = 1.0/(pvalues[testIndex][stepIndex].size()*1.0);

        int npvals = pvalues[testIndex][stepIndex].size();
        for(int i = 0; i < npvals; i++){
            qq.push_back(-std::log10(1-(i * step)));
            p.push_back(-std::log10(pvalues[testIndex][stepIndex][npvals-(i+1)]));
        }

        ui->simplot_plot2->addGraph();
        ui->simplot_plot2->graph()->setData(qq, p);
        QPen scatterPen = QPen(colours[testIndex], 2);

        if(focusGraph >= 0 && focusGraph != testIndex)
            scatterPen = QPen(QColor(colours[testIndex].red(),
                                     colours[testIndex].green(),
                                     colours[testIndex].blue(), 40), 2);

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
    horizontal->setLayer(ui->simplot_plot2->layer("overlay"));
    horizontal->setPen(hpen);
    horizontal->start->setCoords(-1, -std::log10(alpha));
    horizontal->end->setCoords(1e4,-std::log10(alpha));
    alphaLine = horizontal;

    ui->simplot_plot2->replot();
}

void SimPlotWindow::on_simplot_alphaDial_valueChanged(int value){
    alpha = std::pow(10, -value/1000.0);
    ui->simplot_alphaTxt->setText(QString::number(alpha));
    buildPlot();
    updateAlphaLine();
    updatePowerValues(stepIndexForPlot2);
}

void SimPlotWindow::on_simplot_alphaTxt_textChanged(const QString &arg1){

    bool ok = false;
    double x = arg1.toDouble(&ok);
    if(ok && x >= 0 && x <=1){
        alpha = x;
        buildPlot();
        updateAlphaLine();
        updatePowerValues(stepIndexForPlot2);
    }
}


