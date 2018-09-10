#include "SimPlotWindow.h"
#include "ui_simplotwindow.h"

void SimPlotWindow::buildPowerPlot(){

    ui->simplot_power->clearPlottables();
    ui->simplot_power->addGraph();

    QVector<double> index;
    for(int i = 1; i <= nsteps; i++)
        index.push_back(i);


    for(int i = 0; i < testsPerStep; i++){
        QVector<double> power;
        for(int j = 0; j < nsteps; j++)
            power.push_back(calculatePower(j*testsPerStep + i, alpha));

        ui->simplot_power->graph()->setData(index, power);
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
    ui->simplot_power->xAxis->setRange(0.95, ui->simplot_power->xAxis->range().upper + 0.05);

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
