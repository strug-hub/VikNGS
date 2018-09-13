#include "PlotWindow.h"
#include "ui_plotwindow.h"

void PlotWindow::mouseMoveWindow(QMouseEvent *event){

    if(!ui->plot_chrPlt->underMouse()){
        updateVariantInfo(focusedVar);

        int datapoints = ui->plot_chrPlt->graph()->dataCount();

        if(datapoints > 1 || (focusedVar < 0 && datapoints > 0))
            updateVariantHighlightLayer(-1);
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
    double min = chromosomes[focusedChr].getMinPos();
    double buffer = (max-min) * 0.05;
    min = min - buffer;
    max = max + buffer;

    if(adjusted.lower < min && adjusted.upper > max){
        adjusted.lower = min;
        adjusted.upper = max;
    }

    if(adjusted.lower < min){
        adjusted.lower = min;
        if(!isZoom)
            adjusted.upper = oldSize + min;
    }

    if(adjusted.upper > max){
        adjusted.upper = max;
        if(!isZoom)
            adjusted.lower = max - oldSize;
    }

    moveRectangle(zoomRect, focusedChr,
                  std::max(adjusted.lower,min),
                  std::min(max, adjusted.upper));
    ui->plot_genomePlt->layer("overlay")->replot();

    ui->plot_chrPlt->xAxis->setRange(adjusted);
}

void PlotWindow::mouseClickChromosome(QMouseEvent *event){

    double coordx = ui->plot_chrPlt->xAxis->pixelToCoord(event->pos().x());
    double coordy = ui->plot_chrPlt->yAxis->pixelToCoord(event->pos().y());

    if(coordy > 0 && coordx > 0){

        double x = event->pos().x();
        double y = event->pos().y();
        int closest = findClosestVariant(x, y, 100);

        if(closest >= 0){

            if (event->button() == Qt::RightButton)
                focusedVar = -1;
            else
                focusedVar = closest;

            updateVariantHighlightLayer(closest);
            updateVariantInfo(closest);
        }
    }
}

int PlotWindow::findClosestVariant(double x, double y, double maxDist){

    double dist;
    double minDist = 9999999999;
    int minIndex = -1;

    for(int i = 0; i < chromosomes[focusedChr].size(); i++){
        Variant* v = chromosomes[focusedChr].getVariant(i);

        dist = pow(ui->plot_chrPlt->yAxis->coordToPixel(chromosomes[focusedChr].getPval(i)) - y, 2) +
                pow(ui->plot_chrPlt->xAxis->coordToPixel(v->getPosition()) - x, 2);
        if(dist < minDist){
            minDist = dist;
            minIndex = i;
        }
    }

    if(minDist < maxDist)
        return minIndex;

    return -1;
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
    int closest = -1;

    if(coordy > 0 && coordx > 0){
        double x = event->pos().x();
        double y = event->pos().y();

        closest = findClosestVariant(x, y, 100);

        if(closest >= 0)
            updateVariantInfo(closest);
    }

    updateVariantHighlightLayer(closest);
}
