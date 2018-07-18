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

    void on_table_variantTbl_cellClicked(int row, int column);

private:
    Ui::TableDisplayWindow *ui;
    std::vector<Variant>* variants;
    SimulationRequest* request;
    VectorXd y;
    VectorXd g;
    std::vector<CheckTree*> variantCheckTree;

};

#endif // TABLEDISPLAYWINDOW_H
