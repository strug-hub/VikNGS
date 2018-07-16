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

    void initialize(std::vector<std::vector<Variant>>& variants, SimulationRequest& req, QString title);
    void getPvalues(std::vector<std::vector<Variant>>& variants);

    std::vector<std::vector<Variant>> filterCollapsed(std::vector<std::vector<Variant>> &variants, int k);
    void buildPlot();
    void buildLegend();

    void updateSampleSize(int index);

    void mouseMovePlot1(QMouseEvent *event);
    void mouseMovePlot2(QMouseEvent *event);
    void mouseClickPlot1(QMouseEvent *event);
    void mouseClickPlot2(QMouseEvent *event);

private slots:

    void updateGenotypeTable(int index);
    void on_simplot_alphaDial_valueChanged(int value);
    void on_simplot_alphaTxt_textChanged(const QString &arg1);

private:
    Ui::SimPlotWindow *ui;
    QColor grey1 = QColor::fromRgb(190, 190, 190);
    QColor grey2 = QColor::fromRgb(169, 169, 169);
    QColor highlight = QColor::fromRgb(255, 127, 80);
    QColor focus = QColor::fromRgb(102, 204, 204);

    QFont axisFont = QFont("sans", 10, QFont::Bold);

    std::vector<std::vector<Variant>> variants;
    std::vector<std::vector<std::vector<double>>> pvalues;

    SimulationRequest request;
    double alpha;
    int ntests;
    int nsteps;

    QString yAxisLabel;
    QString xAxisLabel;
    QVector<QColor> colours;
    QColor redLine = QColor(210, 80, 80, 125);
    QVector<QString> testTypes;
    int stepIndexForPlot2;
    int testIndexForPlot2;
    QVector<double> calculatePower(int testIndex, double alpha);
    int findClosestPoint(QCustomPlot *plot, QMouseEvent *event, bool getGraphIndex=false);
    void buildPlot2(int stepIndex, int focusGraph = -1);

    QCPItemLine *alphaLine;
    void updateAlphaLine();


};

#endif // SIMPLOTWINDOW_H
