#include "PlotWindow.h"
#include "ui_plotwindow.h"


void PlotWindow::mouseMoveWindow(QMouseEvent *event){

    if(!ui->plot_chrPlt->underMouse()){
        updateVariantInfo(focusedVar);

        int datapoints = ui->plot_chrPlt->graph()->dataCount();

        if(datapoints > 1 || (!focusedVar.isValid() && datapoints > 0))
            updateVariantHighlightLayer(nullVariant);
    }

    if(!ui->plot_genomePlt->underMouse()){
        if(highlightChr != ""){

            moveRectangle(focusRect, focusedChr);
            highlightChr = "";
            ui->plot_genomePlt->layer("overlay")->replot();

        }
    }
}

void PlotWindow::mouseClickGenome(QMouseEvent *event){

    QString chrName = getChromUnderCursor(event);

    if(chrName == focusedChr || chrName == "")
        return;

    if(focusedChr != "")
        resetColor(focusedChr);

    focusedChr = chrName;

    getGraphByName(ui->plot_genomePlt, chrName)->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                focus, Qt::white, 2));

    moveRectangle(focusRect, chrName);

    ui->plot_genomePlt->replot();

    if(chrName != "")
        buildChromosomePlot(chrName);

}

void PlotWindow::updateChromosomeRange(const QCPRange &newRange){

    QCPRange r = ui->plot_chrPlt->xAxis->range();
    double oldSize = r.upper - r.lower;
    double newSize = newRange.upper - newRange.lower;

    bool isZoom = std::abs(newSize - oldSize) > 0.0001;

    QCPRange adjusted = QCPRange();
    adjusted.upper = newRange.upper;
    adjusted.lower = newRange.lower;

    double max = chromosomes[focusedChr].getMaxPos();
    max = max + max * 0.01;

    if(newRange.lower < 0){
        adjusted.lower = 0;
        if(!isZoom)
            adjusted.upper = std::min(max, oldSize);
    }
    if(newRange.upper > max){
          adjusted.upper = max;
          if(!isZoom)
              adjusted.lower = std::max(0.0, max-oldSize);
    }

    moveRectangle(zoomRect, focusedChr, adjusted.lower, adjusted.upper);
    ui->plot_genomePlt->layer("overlay")->replot();

    ui->plot_chrPlt->xAxis->setRange(adjusted);
}

void PlotWindow::mouseClickChromosome(QMouseEvent *event){

    double coordx = ui->plot_chrPlt->xAxis->pixelToCoord(event->pos().x());
    double coordy = ui->plot_chrPlt->yAxis->pixelToCoord(event->pos().y());

    if(coordy > 0 && coordx > 0){

        double x = event->pos().x();
        double y = event->pos().y();
        Variant closest = findClosestVariant(x, y, 100);

        if(closest.isValid()){

            if (event->button() == Qt::RightButton)
                focusedVar = nullVariant;
            else
                focusedVar = closest;

            updateVariantHighlightLayer(closest);
            updateVariantInfo(closest);
        }
    }
}

Variant PlotWindow::findClosestVariant(double x, double y, double maxDist){

    double dist;
    double minDist = 9999999999;
    int minIndex;

    for(int i = 0; i < chromosomes[focusedChr].size(); i++){
        Variant v = chromosomes[focusedChr].getVariant(i);

        dist = pow(ui->plot_chrPlt->yAxis->coordToPixel(-log10(v.getPval(0))) - y, 2) +
                pow(ui->plot_chrPlt->xAxis->coordToPixel(v.pos) - x, 2);
        if(dist < minDist){
            minDist = dist;
            minIndex = i;
        }
    }

    if(minDist < maxDist)
        return chromosomes[focusedChr].getVariant(minIndex);

    return nullVariant;
}

void PlotWindow::mouseMoveGenome(QMouseEvent *event){

    QString chrName = getChromUnderCursor(event);
    if(highlightChr == chrName)
        return;
    if(chrName != "")
        moveRectangle(focusRect, chrName);
    else
        moveRectangle(focusRect, focusedChr);

    highlightChr = chrName;
    ui->plot_genomePlt->layer("overlay")->replot();
}

void PlotWindow::mouseMoveChromosome(QMouseEvent *event){

    double coordx = ui->plot_chrPlt->xAxis->pixelToCoord(event->pos().x());
    double coordy = ui->plot_chrPlt->yAxis->pixelToCoord(event->pos().y());

    updateVariantInfo(focusedVar);
    Variant closest = nullVariant;

    if(coordy > 0 && coordx > 0){
        double x = event->pos().x();
        double y = event->pos().y();

        closest = findClosestVariant(x, y, 100);

        if(closest.isValid())
            updateVariantInfo(closest);
    }

    updateVariantHighlightLayer(closest);

     return;
}
