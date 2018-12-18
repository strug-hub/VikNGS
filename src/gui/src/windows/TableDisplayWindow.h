#ifndef TABLEDISPLAYWINDOW_H
#define TABLEDISPLAYWINDOW_H

#include <QTableWidgetItem>
#include "../Variant.h"
#include "../vikNGS.h"
#include "../simulation/Simulation.h"

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

    void initialize(QString title, Data* data, std::vector<int> testIndexes);

private slots:
    void fillGenotypeTable(int variantIndex);
    void buildVariantTable();
    void on_table_variantTbl_cellClicked(int row, int column);

private:
    Ui::TableDisplayWindow *ui;
    int selectedVariantIndex = 0;

    std::vector<TestSettings>* tests;
    std::vector<int> testIndexes;
    std::vector<VariantSet>* variants;
    int nsamples;
    VectorXd y;
    VectorXi g;
    Family family;

    //variants table
    void addInfo(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table);
    void addPvals(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table);
    void addMafs(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table);
    void addMafsCaseControl(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, bool useCases=true);
    double calculateMaf(VectorXd* gt, bool caseOnly=false, bool controlOnly=false);

    //genotype table
    void addSampleInfo(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table);
    void addSampleGenotypes(int nrow, QStringList &titles, QVector<QVector<QTableWidgetItem*>>& table, Variant* variant);

};

#endif // TABLEDISPLAYWINDOW_H
