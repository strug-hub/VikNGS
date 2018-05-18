#ifndef QZOOMBAR_H
#define QZOOMBAR_H

#include <QWidget>
#include <QPainter>
#include <QMouseEvent>

class QZoomBar : public QWidget
{
    Q_OBJECT

public:
    QZoomBar(QWidget *parent=0);
    void updateMin(double min){ minValue = min; minZoomValue = min;}
    void updateMax(double max){ maxValue = max; maxZoomValue = max;}
    void updateZoomMin(double min){ minZoomValue = min; }
    void updateZoomMax(double max){ maxZoomValue = max; }

signals:
    void rangeChanged(double min, double max);

private:
    double minValue = -1;
    double maxValue = -1;
    int lastMouseX = -1;
    double minZoomValue = -1;
    double maxZoomValue = -1;
    QColor background = QColor::fromRgb(201, 216, 216);
    QColor bar = QColor::fromRgb(102, 204, 204);
    bool mouseDown = false;

    double range(){ return maxValue - minValue; }
    double zoomRange(){ return minZoomValue - maxZoomValue; }

protected:
    virtual void paintEvent(QPaintEvent *event) override;
    virtual void mouseMoveEvent(QMouseEvent *event) override;
    virtual void mouseReleaseEvent(QMouseEvent *event) override { mouseDown = false; }
    virtual void mousePressEvent(QMouseEvent *event) override;

};

#endif // QZOOMBAR_H
