#include "plotwindow.h"
#include "ui_plotwindow.h"


void PlotWindow::mouseMoveWindow(QMouseEvent *event){

    if(!ui->plot_chrPlt->underMouse()){
        updateVariantInfo(focusedVar);

        int datapoints = ui->plot_chrPlt->graph()->dataCount();

        if(datapoints > 1 || (!focusedVar.isValid() && datapoints > 0))
            updateVariantHighlightLayer(nullVariant);

    }

    if(!ui->plot_genomePlt->underMouse()){

        for(int i = 0; i< ui->plot_genomePlt->graphCount(); i++){
            QColor c = ui->plot_genomePlt->graph(i)->scatterStyle().pen().color();
            if( c== highlight)
                mouseMoveGenome(event);
        }

    }

}

void PlotWindow::mouseClickGenome(QMouseEvent *event){

    QString chrName = getChromUnderCursor(event);

    if(chrName != ""){
        focusedVar = nullVariant;
        buildChromosomePlot(chrName, -1, -1);
        mouseMoveGenome(event);
    }   
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

        dist = pow(ui->plot_chrPlt->yAxis->coordToPixel(-log10(v.pvalue)) - y, 2) +
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

    for(int i = 0; i < chrNames.size(); i++){
        getGraphByName(ui->plot_genomePlt, chrNames[i])->setScatterStyle(
                    QCPScatterStyle(QCPScatterStyle::ssDisc,
                                    chromosomes[chrNames[i]].getColour(), Qt::white, 2));
    }

    QString chrName = getChromUnderCursor(event);

    if(chrName != "")
        getGraphByName(ui->plot_genomePlt, chrName)->setScatterStyle(
                    QCPScatterStyle(QCPScatterStyle::ssDisc,
                                    highlight, Qt::white, 2));

    getGraphByName(ui->plot_genomePlt, focusedChr)->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                focus, Qt::white, 2));

    ui->plot_genomePlt->replot();
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

void PlotWindow::mouseScrollChromosome(QWheelEvent *event){
    double scrollSpeed = 0.1;

    double x = ui->plot_chrPlt->xAxis->pixelToCoord(event->pos().x());
    double nsnp = getGraphByName(ui->plot_chrPlt, chrGraphName)->dataCount();

    double minX = ui->plot_chrPlt->xAxis->range().lower;
    double maxX = ui->plot_chrPlt->xAxis->range().upper;
    double range = maxX - minX;

    double fraction = (x-minX)/range;

    if(event->delta() > 0){
        if(nsnp >= 0 && nsnp <=100)
            return;
    double newMin = minX + scrollSpeed*fraction*range;
    double newMax = maxX - scrollSpeed*(1-fraction)*range;

    buildChromosomePlot(focusedChr, newMin, newMax);
    }
    else if(event->delta() < 0){
        double newMin = minX - scrollSpeed*range;
        double newMax = maxX + scrollSpeed*range;

        buildChromosomePlot(focusedChr, newMin, newMax);
    }
}
