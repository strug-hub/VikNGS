#ifndef SIMPLOTWINDOW_H
#define SIMPLOTWINDOW_H

#include <QWidget>
#include "../src/Variant.h"
#include "../widgets/qcustomplot.h"
#include "../simulation/Simulation.h"
#include "../widgets/qcpdocumentobject.h"

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

    void initialize(Data& results, SimulationRequest& req, QString title);

    void buildPowerPlot();
    void buildLegend();

    void updateSampleTable(int index);

    void mouseMovePlot1(QMouseEvent *event);
    void mouseMovePlot2(QMouseEvent *event);
    void mouseClickPlot1(QMouseEvent *event);
    void mouseClickPlot2(QMouseEvent *event);

private slots:

    void on_simplot_alphaDial_valueChanged(int value);
    void on_simplot_alphaTxt_textChanged(const QString &arg1);
    void on_pushButton_clicked(bool checked);

    void on_simplot_pdfBtn_pressed();
    void saveAsPdf(QCustomPlot* plot, QString fileName, int stepIndex=-1);

private:
    Ui::SimPlotWindow *ui;
    QColor grey1 = QColor::fromRgb(190, 190, 190);
    QColor grey2 = QColor::fromRgb(169, 169, 169);
    QColor highlight = QColor::fromRgb(255, 127, 80);
    QColor focus = QColor::fromRgb(102, 204, 204);
    QFont axisFont = QFont("sans", 10, QFont::Bold);
    QColor redLine = QColor(210, 80, 80, 125);

    QString yAxisLabel;
    QString xAxisLabel;

    SimulationRequest request;
    Data result;

    void getPvalues(std::vector<VariantSet>& variants);
    std::vector<std::vector<double>> pvalues;
    std::vector<Test> tests;

    int nsteps;
    int testsPerStep;
    inline int ntests() { return nsteps * testsPerStep; }

    QVector<QColor> colours;
    int powerIndex = 0;
    int qqIndex = 0;
    double calculatePower(int testIndex, double alpha);

    int findClosestPoint(QCustomPlot *plot, QMouseEvent *event, bool getGraphIndex=false);
    void buildQQPlot(int stepIndex, int focusGraph = -1);
    void updatePowerValues(int index);

    QCPItemLine *alphaLine;
    double alpha;
    void updateAlphaLine();

    QString lastSaveDir = ".";

};

#endif // SIMPLOTWINDOW_H
