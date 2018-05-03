#include "qzoombar.h"
#include <algorithm>

QZoomBar::QZoomBar(QWidget *parent) : QWidget(parent) {
    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);
}

void QZoomBar::paintEvent(QPaintEvent * event) {
    // call the base class paint event method
    // this will draw the base class content
    QWidget::paintEvent(event);
    // draw a blue border inside the button
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen(QPen(background, 2));

    const int radius = height()/2.0;

    QRect rectangleFull = QRect(0, 0, width(), height());
    QPainterPath bgPath;
    bgPath.addRoundedRect(rectangleFull, radius, radius);
    painter.setPen(QPen(background, 2));
    painter.fillPath(bgPath, QBrush(background, Qt::SolidPattern));
    painter.drawPath(bgPath);

    const double propor1 = std::max(0.0, (minZoomValue - minValue)/range());
    const double propor2 = std::max(0.0, (maxValue - maxZoomValue)/range());

    const int x = width()*propor1;
    const int w = width() - width()*propor1 - width()*propor2;

    QRect rectangleZoom = QRect(x, 0, w, height());
    QPainterPath barPath;
    barPath.addRoundedRect(rectangleZoom, radius, radius);
    painter.setPen(QPen(bar, 2));
    painter.fillPath(barPath, QBrush(bar, Qt::Dense2Pattern));
    painter.drawPath(barPath);

    //painter.drawRoundedRect(rectangleZoom,2,2);
    //painter.fillRect(rectangleZoom, QBrush(bar, Qt::Dense2Pattern));

    this->update();
}


void QZoomBar::mouseMoveEvent(QMouseEvent * event) {

    int delta = lastMouseX - event->pos().x();

    double scale = range()/width();

    if(minZoomValue <= minValue && maxZoomValue >= maxValue)
        return;

    if(delta > 0){
        double oldX = minZoomValue;
        minZoomValue = std::max(minValue, minZoomValue - delta*scale);
        maxZoomValue -= oldX - minZoomValue;
    }
    if(delta < 0){
        double oldX = maxZoomValue;
        maxZoomValue = std::min(maxValue, maxZoomValue - delta*scale);
        double shift = maxZoomValue - oldX;
        minZoomValue += shift;
    }

    emit rangeChanged(minZoomValue, maxZoomValue);

    lastMouseX = event->pos().x();

}


void QZoomBar::mousePressEvent(QMouseEvent *event) {
    int mouseX = event->pos().x();

    const double propor1 = std::max(0.0, (minZoomValue - minValue)/range());
    const double propor2 = std::max(0.0, (maxValue - maxZoomValue)/range());

    const int x1 = width()*propor1;
    const int x2 = width() - width()*propor2;

    if(mouseX > x1 && mouseX < x2)
        lastMouseX = event->pos().x();
    else{
        lastMouseX = ((maxZoomValue - minZoomValue)/2 + minZoomValue) / (range()/width());
        mouseMoveEvent(event);
    }

    mouseDown = true;
}
