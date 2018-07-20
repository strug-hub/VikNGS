#ifndef TABLEDISPLAYWINDOW_H
#define TABLEDISPLAYWINDOW_H

#include <QWidget>
#include "./checktree.h"
#include "../src/Variant.h"
#include "../simulation/simulation.h"

namespace Ui {
class TableDisplayWindow;
}

class TableDisplayWindow : public QWidget
{
    Q_OBJECT

public:
    explicit TableDisplayWindow(QWidget *parent = 0);
    ~TableDisplayWindow();

public slots:

    void initialize(QString title, std::vector<Variant>& variants, SimulationRequest& request, int index);

private slots:
    void fillGenotypeTable(int variantIndex);
    void buildVariantTable();
    void drawVariantCheckTree();
    void displayCalls(int index);

    void on_table_variantTbl_cellClicked(int row, int column);
    void on_table_genoTbl_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn);

private:
    Ui::TableDisplayWindow *ui;
    std::vector<Variant>* variants;
    int selectedVariantIndex = 0;

    SimulationRequest* request;
    VectorXd y;
    VectorXd g;
    std::vector<CheckTree*> variantCheckTree;

    QColor colourA = QColor(0,151,45);
    QColor colourT = QColor(255,0,43);
    QColor colourC = QColor(0,37,248);
    QColor colourG = QColor(222,111,47);

    QString basePair(QString b){
        QString html = "<font color=\"";

        if(b == "A")
            html.append(colourA.name());
        else if (b == "C")
            html.append(colourC.name());
        else if (b == "G")
            html.append(colourG.name());
        else if (b == "T")
             html.append(colourT.name());

        html.append("\">").append(b).append("</font>");
        return html;
    }

};

#endif // TABLEDISPLAYWINDOW_H
