#include "SimPlotWindow.h"
#include "ui_simplotwindow.h"
#include "VikngsUiConstants.h"

void SimPlotWindow::buildPowerPlot(){

    ui->simplot_power->clearPlottables();
    ui->simplot_power->addGraph();

    QVector<double> size;
    for(int i = 1; i <= nsteps; i++){
        int n = 0;
        for(SimulationRequestGroup srg : request.groups)
            n+=srg.getSampleSize(i);

        size.push_back(n);
    }

    for(int i = 0; i < testsPerStep; i++){
        QVector<double> power;
        for(int j = 0; j < nsteps; j++)
            power.push_back(calculatePower(j*testsPerStep + i, alpha));

        ui->simplot_power->graph()->setData(size, power);
        ui->simplot_power->graph()->setPen(QPen(colours[i], 2));
        ui->simplot_power->graph()->setScatterStyle(QCPScatterStyle::ssDisc);
        //last layer is highlight layer
        ui->simplot_power->addGraph();
    }

    ui->simplot_power->rescaleAxes();
    ui->simplot_power->xAxis->setLabel(xAxisLabel);
    ui->simplot_power->yAxis->setLabel(yAxisLabel);
    ui->simplot_power->xAxis->setLabelFont(axisFont);
    ui->simplot_power->yAxis->setLabelFont(axisFont);
    ui->simplot_power->yAxis->setRange(0, 1.05);

    float minXSize = size[0];
    float maxXSize = size.back();
    float xRange = maxXSize - minXSize;
    ui->simplot_power->xAxis->setRange(minXSize - 0.025*xRange, maxXSize + 0.025*xRange);

    ui->simplot_power->replot();
}

void SimPlotWindow::buildLegend(){

    if(testsPerStep > 0){
        QString t = QString::fromStdString(result.tests[0].toString());
        ui->simplot_legend1Lbl->setText(t);
        QPalette pal = palette();
        pal.setColor(QPalette::Background, colours[0]);
        ui->simplot_legend1Sqr->setAutoFillBackground(true);
        ui->simplot_legend1Sqr->setPalette(pal);
    }
    else{
        ui->simplot_legend1Lbl->hide();
        ui->simplot_legend1Sqr->hide();
    }

    if(testsPerStep > 1){
        QString t = QString::fromStdString(result.tests[1].toString());
        ui->simplot_legend2Lbl->setText(t);
        QPalette pal = palette();
        pal.setColor(QPalette::Background, colours[1]);
        ui->simplot_legend2Sqr->setAutoFillBackground(true);
        ui->simplot_legend2Sqr->setPalette(pal);
    }
    else{
        ui->simplot_legend2Lbl->hide();
        ui->simplot_legend2Sqr->hide();
    }

    if(testsPerStep > 2){
        QString t = QString::fromStdString(result.tests[2].toString());
        ui->simplot_legend3Lbl->setText(t);
        QPalette pal = palette();
        pal.setColor(QPalette::Background, colours[2]);
        ui->simplot_legend3Sqr->setAutoFillBackground(true);
        ui->simplot_legend3Sqr->setPalette(pal);
    }
    else{
        ui->simplot_legend3Lbl->hide();
        ui->simplot_legend3Sqr->hide();
    }

    if(testsPerStep > 3){
        QString t = QString::fromStdString(result.tests[3].toString());
        ui->simplot_legend4Lbl->setText(t);
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

    if(powerIndex >= 0){
        alphaLine->start->setCoords(-1, -std::log10(alpha));
        alphaLine->end->setCoords(1e4,-std::log10(alpha));

        ui->simplot_plot2->layer("overlay")->replot();
    }
}

void SimPlotWindow::buildQQPlot(int stepNumber, int focusGraph){

    if(stepNumber < 0)
        return;

    if(focusGraph < 0)
        focusGraph = qqIndex;

    ui->simplot_plot2->clearItems();
    ui->simplot_plot2->clearPlottables();

    double max = 0;
    int start = stepNumber * testsPerStep;

    for(int i = start; i < start + testsPerStep; i++){

        QVector<double> qq;
        QVector<double> p;

        int npvals = static_cast<int>(pvalues[i].size());
        double expectedIncrement = 1.0/(npvals*1.0);

        for(int j = 0; j < npvals; j++){
            qq.push_back(-std::log10(1-(j * expectedIncrement)));
            p.push_back(-std::log10(pvalues[i][npvals - j - 1]));
        }

        ui->simplot_plot2->addGraph();
        ui->simplot_plot2->graph()->setData(qq, p);
        int step = i - start;
        QPen scatterPen = QPen(colours[step], 2);

        if(focusGraph >= 0 && focusGraph != step)
            scatterPen = QPen(QColor(colours[step].red(),
                                     colours[step].green(),
                                     colours[step].blue(), 40), 2);

        ui->simplot_plot2->graph()->setPen(scatterPen);
        ui->simplot_plot2->graph()->setScatterStyle(QCPScatterStyle::ssDisc);
        ui->simplot_plot2->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);

        double testMax = std::max(p.last(), qq.last());
        if(testMax > max)
            max = testMax;
    }

    if(focusGraph >= 0){
        if(ui->simplot_plot2->layer(TOP_LAYER) == ui->simplot_plot2->layer(NO_LAYER))
            ui->simplot_plot2->addLayer(TOP_LAYER, ui->simplot_plot2->layer(MAIN_LAYER),
                                        QCustomPlot::LayerInsertMode::limAbove);
        ui->simplot_plot2->graph(focusGraph)->setLayer(TOP_LAYER);
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
    horizontal->setLayer(ui->simplot_plot2->layer(OVERLAY_LAYER));
    horizontal->setPen(hpen);
    horizontal->start->setCoords(-1, -std::log10(alpha));
    horizontal->end->setCoords(1e4,-std::log10(alpha));
    alphaLine = horizontal;

    ui->simplot_plot2->replot();
    buildHistPlot(stepNumber, focusGraph);
}


QPair<QVector<double>, QVector<double>> histBin(std::vector<double>& values, double nbins, double min, double max){

    QVector<double> bin;
    QVector<double> height;

    double binSize = (max - min)/nbins;

    double width = min;
    while(width < max){
        width += binSize;
        bin.push_back(width);
        height.push_back(0);
    }

    for(int i = 0; i < values.size(); i++){
        for(int j = 0; j< bin.size(); j++){
            if(values[i] < bin[j]){
                height[j]++;
                j = bin.size() + 1;
            }
        }
    }

    for(int i = 0; i < height.size(); i++){
        height[i] /= values.size();
        bin[i] -= binSize/2.0;
    }

    QPair<QVector<double>, QVector<double>> ret;
    ret.first = bin;
    ret.second = height;
    return ret;
}

void SimPlotWindow::buildHistPlot(int stepNumber, int focusGraph){

    int nbins = ui->simplot_histBinBox->value();

    QCustomPlot* plt = ui->simplot_histPlot;

    plt->clearItems();
    plt->clearPlottables();
    plt->setBackground(Qt::transparent);
    plt->setAttribute(Qt::WA_OpaquePaintEvent, false);

    if(plt->layer(TOP_LAYER) == plt->layer(NO_LAYER)){
        plt->addLayer(TOP_LAYER, plt->layer(MAIN_LAYER),
                                    QCustomPlot::LayerInsertMode::limAbove);
    }

    double max = 0;
    int start = stepNumber * testsPerStep;

    for(int i = start; i < start + testsPerStep; i++){

        QPair<QVector<double>, QVector<double>> binHeights;
        binHeights = histBin(pvalues[i], nbins, 0, 1);
        plt->addGraph();
        QCPBars *bars = new QCPBars(plt->xAxis, plt->yAxis);
        bars->setWidth(1.0/nbins);
        bars->setPen(NO_PEN);
        bars->setData(binHeights.first, binHeights.second);

        int step = i - start;
        if(focusGraph >= 0 && focusGraph == step){
            plt->graph()->setLayer(TOP_LAYER);
            bars->setLayer(TOP_LAYER);
            bars->setBrush(brush(colours[step]));
            bars->setPen(pen(BLACK, 0.4));
        }
        else
            bars->setBrush(brush(colours[step], 0.4));

        double maxHeight = *std::max_element(binHeights.second.constBegin(), binHeights.second.constEnd());
        if(maxHeight > max)
            max = maxHeight;
    }


    QPen black_25 = pen(BLACK, 0.25);
    QPen black_50 = pen(BLACK, 0.5);

    plt->xAxis->setBasePen(black_25);
    plt->yAxis->setBasePen(black_25);
    plt->xAxis->setTickLabelColor(black_50.color());
    plt->yAxis->setTickLabelColor(black_50.color());
    plt->xAxis->setTickPen(black_25);
    plt->yAxis->setTickPen(black_25);
    plt->xAxis->setSubTickPen(black_25);
    plt->yAxis->setSubTickPen(black_25);

    plt->xAxis->grid()->setPen(NO_PEN);
    plt->yAxis->grid()->setPen(NO_PEN);

    QFont pfont = plt->xAxis->tickLabelFont();
    pfont.setPointSize(6);
    plt->xAxis->setTickLabelFont(pfont);
    plt->yAxis->setTickLabelFont(pfont);

    plt->yAxis->setRange(0, max * 1.05);
    plt->xAxis->setRange(0, 1);

    plt->replot();
}
