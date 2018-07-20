#ifndef TABLEDISPLAYWINDOW_H
#define TABLEDISPLAYWINDOW_H

#include <QTableWidgetItem>
#include "./checktree.h"
#include "../src/Variant.h"
#include "../src/RVS.h"
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

    void initialize(QString title, std::vector<Variant>& variants, SimulationRequest& simRequest, int index);
    void initialize(QString title, Result* result);

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

    SimulationRequest* simRequest;

    bool simulation = true;

    VectorXd y;
    VectorXd g;
    std::vector<CheckTree*> variantCheckTree;

    QColor colourA = QColor(0,151,45);
    QColor colourT = QColor(255,0,43);
    QColor colourC = QColor(0,37,248);
    QColor colourG = QColor(222,111,47);

    void addInfo(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table);
    void addPvals(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table);
    void addMafs(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, QString gt="expected", QString caseControl="no");
    void addGtFreq(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, QString gt="expected", QString caseControl="no");
    void addSampleInfo(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, Variant& variant);
    void addSampleGenotypes(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, Variant& variant);
    void addVCFData(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, Variant& variant);
    void addAlleleProbability(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, Variant& variant);

    double calculateMaf(VectorXd& gt, bool caseOnly=false, bool controlOnly=false);

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
