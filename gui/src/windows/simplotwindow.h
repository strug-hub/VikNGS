#ifndef SIMPLOTWINDOW_H
#define SIMPLOTWINDOW_H

#include <QWidget>
#include "../src/Variant.h"
#include "../widgets/qcustomplot.h"
#include "../simulation/simulation.h"


namespace Ui {
class SimPlotWindow;
}

class SimPlotWindow : public QWidget
{
    Q_OBJECT

public:
    explicit SimPlotWindow(QWidget *parent = 0);
    ~SimPlotWindow();

public slots:

    void initialize(std::vector<std::vector<Variant>> variants, std::vector<SimulationRequest> reqs, QString title);
    void buildPlot();
    void mouseMovePlot1(QMouseEvent *event);
    void mouseClickPlot1(QMouseEvent *event);

private:
    Ui::SimPlotWindow *ui;
    QColor grey1 = QColor::fromRgb(190, 190, 190);
    QColor grey2 = QColor::fromRgb(169, 169, 169);
    QColor highlight = QColor::fromRgb(255, 127, 80);
    QColor focus = QColor::fromRgb(102, 204, 204);

    std::vector<std::vector<Variant>> variants;
    std::vector<SimulationRequest> requests;

    double calculatePower(std::vector<Variant> run, double alpha);
    QVector<double> calculatePower(std::vector<std::vector<Variant>> variants, double alpha);
    int findClosestPoint(QCustomPlot *plot, QMouseEvent *event);
    void buildPlot2(int index);



};

#endif // SIMPLOTWINDOW_H
