#include "PlotWindow.h"
#include "ui_plotwindow.h"
#include <QSharedPointer>

void PlotWindow::updateVariantHighlightLayer(int variantIndex){
    QVector<double> l;
    QVector<double> p;

    if(variantIndex >= 0){
        l.push_back(chromosomes[focusedChr].getPosition(variantIndex));
        p.push_back(chromosomes[focusedChr].getPval(variantIndex));
    }

    if(focusedVar >= 0){
        l.push_back(chromosomes[focusedChr].getPosition(focusedVar));
        p.push_back(chromosomes[focusedChr].getPval(focusedVar));
    }

    ui->plot_chrPlt->graph()->setData(l,p);
    ui->plot_chrPlt->replot();
}

void PlotWindow::buildGenomePlot(){
    QCustomPlot* plt = ui->plot_genomePlt;
    plt->clearPlottables();
    plt->xAxis->grid()->setPen(Qt::NoPen);
    //plt->xAxis->setTickPen(Qt::NoPen);

    bool colour1 = true;    
    double offset = 0;
    QSharedPointer<QCPAxisTicker> ticks = QSharedPointer<QCPAxisTickerText>(new QCPAxisTickerText);

    for(int i = 0; i < chrNames.size(); i++){
        chromosomes[chrNames[i]].setOffset(offset);
        Chromosome chr = chromosomes[chrNames[i]];

        plt->addGraph();
        plt->graph()->setData(chr.getRelativePositions(offset), chr.getPvals(0));
        plt->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);
        plt->graph()->setName(chr.getName());

        plt->graph()->setAdaptiveSampling(false);
        //plt->graph()->
        QColor toUse = grey1;
        if(!colour1)
            toUse=grey2;

        colour1 = !colour1;

        chromosomes[chrNames[i]].setColour(toUse);

        if(chrNames[i] == focusedChr)
            toUse = focus;

        plt->graph()->setScatterStyle(
                    QCPScatterStyle(QCPScatterStyle::ssDisc,
                                    toUse, Qt::white, 2));

        ((QCPAxisTickerText*)ticks.data())->addTick(offset + chr.getSpan()/2.0, chr.getName());
        offset += chr.getSpan();

        if(i < chrNames.size()-1){
            QCPItemLine *divider;
            QPen vpen = QPen(Qt::SolidLine);
            vpen.setColor( QColor::fromRgb(110, 128, 158, 100));
            divider = new QCPItemLine(plt);
            divider->setPen(vpen);
            divider->start->setCoords(offset,0);
            divider->end->setCoords(offset,1e6);
        }
    }

    plt->xAxis->setTicker(ticks);
    plt->xAxis->setTickLabelRotation(-30);

    focusRect = new QCPItemRect(plt);
    focusRect->setLayer("overlay");
    focusRect->topLeft->setType(QCPItemPosition::ptPlotCoords);
    focusRect->topLeft->setAxisRect( plt->axisRect() );
    focusRect->bottomRight->setType(QCPItemPosition::ptPlotCoords);
    focusRect->bottomRight->setAxisRect( plt->axisRect() );
    focusRect->setPen(QPen(focus));
    QColor transFocus = focus;
    transFocus.setAlpha(50);
    focusRect->setBrush(QBrush(transFocus));

    zoomRect = new QCPItemRect(plt);
    zoomRect->setLayer("overlay");
    zoomRect->topLeft->setType(QCPItemPosition::ptPlotCoords);
    zoomRect->topLeft->setAxisRect( plt->axisRect() );
    zoomRect->bottomRight->setType(QCPItemPosition::ptPlotCoords);
    zoomRect->bottomRight->setAxisRect( plt->axisRect() );
    zoomRect->setPen(Qt::NoPen);
    zoomRect->setBrush(QBrush(transFocus));

    plt->yAxis->setLabel("-log10(p)");

    plt->rescaleAxes();
    plt->yAxis->setRangeLower(0);

    double ymax = plt->yAxis->range().upper;
    if(ymax < 7.5){
        plt->yAxis->setRangeUpper(7.5);
        ymax = 7.5;
    }
    plt->yAxis->setRangeUpper(ymax + ymax * 0.05);

    QPen hpen = QPen(Qt::DashDotLine);
    hpen.setColor( QColor::fromRgb(210, 80, 80));
    horizontal = new QCPItemLine(plt);
    horizontal->setPen(hpen);
    horizontal->start->setCoords(0,1.301);
    horizontal->end->setCoords(plt->xAxis->range().upper,1.301);

    QPen bpen = QPen(Qt::DashDotLine);
    bpen.setColor( QColor::fromRgb(250, 145, 145));
    QCPItemLine *bonferroni = new QCPItemLine(plt);
    bonferroni->setPen(bpen);
    bonferroni->start->setCoords(0,7.301);
    bonferroni->end->setCoords(plt->xAxis->range().upper,7.301);

    moveRectangle(focusRect, plt->graph(0)->name());
    moveRectangle(zoomRect, plt->graph(0)->name());

    plt->replot();
}

void PlotWindow::buildChromosomePlot(QString chrName){
    ui->plot_chrPlt->clearPlottables();

    focusedVar = -1;
    focusedChr = chrName;
    Chromosome chr = chromosomes[chrName];

    ui->plot_chrPlt->addGraph();

    QVector<double> pos = chr.getPositions();
    QVector<double> pvals = chr.getPvals(0);

    ui->plot_chrPlt->graph()->setData(pos, pvals);
    ui->plot_chrPlt->graph()->setName(chrGraphName);
    ui->plot_chrPlt->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);

    double pointSize = 1.0/(sqrt(pos.size())/100.0);
    if(pointSize > MAX_POINT_SIZE) pointSize = MAX_POINT_SIZE;
    if(pointSize < MIN_POINT_SIZE) pointSize = MIN_POINT_SIZE;

    ui->plot_chrPlt->graph()->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                focus, Qt::white, pointSize));

    ui->plot_chrPlt->yAxis->setLabel("-log10(p)");
    ui->plot_chrPlt->xAxis->setLabel("Position");

    QCPTextElement *title = new QCPTextElement(ui->plot_chrPlt, "Chromosome " + chrName, QFont("sans", 17, QFont::Bold));
    ui->plot_chrPlt->plotLayout()->removeAt(0);
    ui->plot_chrPlt->plotLayout()->addElement(0,0,title);

    ui->plot_chrPlt->rescaleAxes();
    ui->plot_chrPlt->yAxis->setRangeLower(0);

    double ymax = ui->plot_chrPlt->yAxis->range().upper;
    if(ymax < 7.5){
        ui->plot_chrPlt->yAxis->setRangeUpper(7.5);
        ymax = 7.5;
    }
    ui->plot_chrPlt->yAxis->setRangeUpper(ymax + ymax * 0.05);

    QPen hpen = QPen(Qt::DashDotLine);
    hpen.setColor( QColor::fromRgb(210, 80, 80));
    horizontal = new QCPItemLine(ui->plot_chrPlt);
    horizontal->setPen(hpen);
    horizontal->start->setCoords(0,1.301);
    horizontal->end->setCoords(chr.getMaxPos()*2,1.301);

    QPen bpen = QPen(Qt::DashDotLine);
    bpen.setColor( QColor::fromRgb(250, 145, 145));
    QCPItemLine *bonferroni = new QCPItemLine(ui->plot_chrPlt);
    bonferroni->setPen(bpen);
    bonferroni->start->setCoords(0,7.301);
    bonferroni->end->setCoords(chr.getMaxPos()*2, 7.301);

    ui->plot_chrPlt->addGraph();
    ui->plot_chrPlt->graph()->setName("annotationLayer");
    ui->plot_chrPlt->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);
    ui->plot_chrPlt->graph()->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                highlight, Qt::white, MAX_POINT_SIZE));

    updateVariantHighlightLayer(-1);

    ui->plot_chrPlt->setInteraction(QCP::iRangeDrag, true);
    ui->plot_chrPlt->setInteraction(QCP::iRangeZoom, true);
    ui->plot_chrPlt->axisRect()->setRangeDrag(Qt::Horizontal);
    ui->plot_chrPlt->axisRect()->setRangeZoom(Qt::Horizontal);


     ui->plot_chrPlt->replot();
}

