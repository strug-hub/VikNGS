#ifndef PLOTWINDOW_H
#define PLOTWINDOW_H

#include <QWidget>
#include "../Variant.h"
#include "../widgets/qcustomplot.h"
#include "./TableDisplayWindow.h"
#include "Chromosome.h"
#include "../widgets/qcpdocumentobject.h"

namespace Ui {
class PlotWindow;
}

class PlotWindow : public QWidget
{
    Q_OBJECT

public:
    explicit PlotWindow(QWidget *parent = 0);
    ~PlotWindow();

public slots:

    void initialize(Data& result, QString title);
  //  void initialize(int n, QString title);
    void createChromosomes(std::vector<VariantSet>& variants);

    void updateVariantInfo(int variantIndex);

    void mouseMoveWindow(QMouseEvent* event);
    void mouseMoveGenome(QMouseEvent* event);
    void mouseClickGenome(QMouseEvent *event);
    void updateChromosomeRange(const QCPRange &newRange);
    QString getChromUnderCursor(QMouseEvent *event);
    void resetColor(QString chrName);
    void mouseMoveChromosome(QMouseEvent *event);
    void mouseClickChromosome(QMouseEvent *event);

    void buildGenomePlot();
    void buildChromosomePlot(QString chrName);
    void updateVariantHighlightLayer(int variantIndex);
    void moveRectangle(QCPItemRect *rect, QString chrName, double lower=-1, double upper=-1);

    QCPGraph* getGraphByName(QCustomPlot *plot, QString name);
    int findClosestVariant(double x, double y, double maxDist);

private slots:
    void on_plot_genotypeBtn_pressed();
    void on_plot_pdfBtn_pressed();

private:
    Ui::PlotWindow *ui;  
    TableDisplayWindow *tableView = new TableDisplayWindow();
    const double MAX_POINT_SIZE = 8;
    double MIN_POINT_SIZE = 2;

    QColor grey1 = QColor::fromRgb(190, 190, 190);
    QColor grey2 = QColor::fromRgb(169, 169, 169);
    QColor highlight = QColor::fromRgb(255, 127, 80);
    QColor focus = QColor::fromRgb(102, 204, 204);

    Data result;

    QMap<QString, Chromosome> chromosomes;
    QVector<QString> chrNames;
    QCPItemLine *horizontal;
    QCPItemRect *focusRect;
    QCPItemRect *zoomRect;

    QString highlightChr = "";
    QString focusedChr = "";

    int focusedVar;
    int highlightVar;

    const QString chrGraphName = "chromgraph";

    QString lastSaveDir = ".";
    void saveAsPdf(QString fileName);

    bool eventFilter(QObject *obj, QEvent *event)
    {
      if (event->type() == QEvent::MouseMove)
      {
        QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
        mouseMoveWindow(mouseEvent);
      }
      return false;
    }
};

#endif // PLOTWINDOW_H
