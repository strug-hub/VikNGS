#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "RVS.h"

#include <QMainWindow>
#include <QFileDialog>
#include <QThread>
#include <QtConcurrent/QtConcurrent>
#include <QMessageBox>
#include <QVector>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:

    void on_sim_runBtn_clicked();
    void on_sim_groupAddBtn_clicked();
    void on_sim_groupRemoveBtn_clicked();

    void sim_replot(QVector<double> values, int maxX);
    std::vector<SimulationRequestGroup> constructGroups(int test, int ntest);
    SimulationRequest constructRequest(std::vector<SimulationRequestGroup> groups);
    double calculatePower( QVector<double> pval);

    void on_sim_testBootChk_stateChanged(int arg1);
    void on_sim_testBootChk_toggled(bool checked);
    void on_sim_testRareCastBtn_toggled(bool checked);
    void on_sim_testRareCalphaBtn_toggled(bool checked);

    //--------------

    void on_main_vcfDirBtn_clicked();
    void on_main_sampleDirBtn_clicked();
    void on_main_bedDirBtn_clicked();

    void main_replot(QVector<double> values);
    void on_main_runBtn_clicked();

    void on_main_testRareCastBtn_toggled(bool checked);
    void on_main_testRareCalphaBtn_toggled(bool checked);
    void on_main_testBootChk_stateChanged(int arg1);
    void on_main_testBootChk_toggled(bool checked);


private:
    Ui::MainWindow *ui;


};

#endif // MAINWINDOW_H
