#include "PlotWindow.h"
#include "ui_plotwindow.h"

void PlotWindow::updateVariantHighlightLayer(VariantSet* variants){
    QVector<double> l;
    QVector<double> p;

    if(variants->isValid()){
        std::vector<Variant>* vPointer = variants->getVariants();
        for(size_t i = 0; i < vPointer->size(); i++){
            l.push_back(vPointer->at(i).getPosition());
            p.push_back(-log10(variants->getPval(0)));
        }
    }
    if(focusedVar->isValid()){
        std::vector<Variant>* vPointer = focusedVar->getVariants();
        for(size_t i = 0; i < vPointer->size(); i++){
            l.push_back(vPointer->at(i).getPosition());
            p.push_back(-log10(focusedVar->getPval(0)));
        }
    }

    ui->plot_chrPlt->graph()->setData(l,p);
    ui->plot_chrPlt->replot();
}

void PlotWindow::buildGenomePlot(){
    ui->plot_genomePlt->clearPlottables();

    QVector<double> tickValues;
    QVector<QString> tickLabels;

    bool colour1 = true;    

    double offset = 0;

    for( int i = 0; i < chrNames.size(); i++ )
    {
        chromosomes[chrNames[i]].setOffset(offset);
        Chromosome chr = chromosomes[chrNames[i]];

        ui->plot_genomePlt->addGraph();
        ui->plot_genomePlt->graph()->setData(chr.getPositions(offset), chr.getPvals(0));
        ui->plot_genomePlt->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);
        ui->plot_genomePlt->graph()->setName(chr.getName());

        tickLabels.push_back(chr.getName());
        tickValues.push_back(offset + chr.getSpan()/2);

        QColor toUse = grey1;
        if(!colour1)
            toUse=grey2;

        colour1 = !colour1;

        chromosomes[chrNames[i]].setColour(toUse);

        if(chrNames[i] == focusedChr)
            toUse = focus;

        ui->plot_genomePlt->graph()->setScatterStyle(
                    QCPScatterStyle(QCPScatterStyle::ssDisc,
                                    toUse, Qt::white, 2));
        offset += chr.getSpan();       
    }

    focusRect = new QCPItemRect(ui->plot_genomePlt);
    focusRect->setLayer("overlay");
    focusRect->topLeft->setType(QCPItemPosition::ptPlotCoords);
    focusRect->topLeft->setAxisRect( ui->plot_genomePlt->axisRect() );
    focusRect->bottomRight->setType(QCPItemPosition::ptPlotCoords);
    focusRect->bottomRight->setAxisRect( ui->plot_genomePlt->axisRect() );
    focusRect->setPen(QPen(focus));
    QColor transFocus = focus;
    transFocus.setAlpha(50);
    focusRect->setBrush(QBrush(transFocus));

    zoomRect = new QCPItemRect(ui->plot_genomePlt);
    zoomRect->setLayer("overlay");
    zoomRect->topLeft->setType(QCPItemPosition::ptPlotCoords);
    zoomRect->topLeft->setAxisRect( ui->plot_genomePlt->axisRect() );
    zoomRect->bottomRight->setType(QCPItemPosition::ptPlotCoords);
    zoomRect->bottomRight->setAxisRect( ui->plot_genomePlt->axisRect() );
    zoomRect->setPen(Qt::NoPen);
    zoomRect->setBrush(QBrush(transFocus));

    ui->plot_genomePlt->yAxis->setLabel("-log(p)");

    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(tickValues, tickLabels);
    ui->plot_genomePlt->xAxis->setTicker(textTicker);
    ui->plot_genomePlt->xAxis->setTickLabelRotation(90);
    ui->plot_genomePlt->xAxis->grid()->setPen(Qt::NoPen);

    ui->plot_genomePlt->rescaleAxes();
    ui->plot_genomePlt->yAxis->setRangeLower(0);

    double ymax = ui->plot_genomePlt->yAxis->range().upper;
    if(ymax < 7.5){
        ui->plot_genomePlt->yAxis->setRangeUpper(7.5);
        ymax = 7.5;
    }
    ui->plot_genomePlt->yAxis->setRangeUpper(ymax + ymax * 0.05);

    QPen hpen = QPen(Qt::DashDotLine);
    hpen.setColor( QColor::fromRgb(210, 80, 80));
    horizontal = new QCPItemLine(ui->plot_genomePlt);
    horizontal->setPen(hpen);
    horizontal->start->setCoords(0,1.301);
    horizontal->end->setCoords(ui->plot_genomePlt->xAxis->range().upper,1.301);

    QPen bpen = QPen(Qt::DashDotLine);
    bpen.setColor( QColor::fromRgb(250, 145, 145));
    QCPItemLine *bonferroni = new QCPItemLine(ui->plot_genomePlt);
    bonferroni->setPen(bpen);
    bonferroni->start->setCoords(0,7.301);
    bonferroni->end->setCoords(ui->plot_genomePlt->xAxis->range().upper,7.301);

    moveRectangle(focusRect, ui->plot_genomePlt->graph(0)->name());
    moveRectangle(zoomRect, ui->plot_genomePlt->graph(0)->name());

    ui->plot_genomePlt->replot();
}

void PlotWindow::buildChromosomePlot(QString chrName){
    ui->plot_chrPlt->clearPlottables();

    focusedVar = &nullVariant;
    focusedChr = chrName;
    Chromosome chr = chromosomes[chrName];

    ui->plot_chrPlt->addGraph();

    QVector<double> pos = chr.getPositions();
    QVector<double> pvals = chr.getPvals(0);

    QVector<double> pos_range;
    QVector<double> pvals_range;

    for(int i =0; i<pos.size(); i++){
        pos_range.push_back(pos[i]);
        pvals_range.push_back(pvals[i]);
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
                                highlight, Qt::white, 6));

    updateVariantHighlightLayer(&nullVariant);

    ui->plot_chrPlt->setInteraction(QCP::iRangeDrag, true);
    ui->plot_chrPlt->setInteraction(QCP::iRangeZoom, true);
    ui->plot_chrPlt->axisRect()->setRangeDrag(Qt::Horizontal);
    ui->plot_chrPlt->axisRect()->setRangeZoom(Qt::Horizontal);


     ui->plot_chrPlt->replot();
}

