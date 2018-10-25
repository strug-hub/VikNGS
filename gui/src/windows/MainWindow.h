#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "../src/vikNGS.h"
#include "../src/Variant.h"
#include "../src/Log.h"
#include "../simulation/Simulation.h"
#include "../AsyncJob.h"

#include <QMainWindow>
#include <QFileDialog>
#include <QThread>
#include <QtConcurrent/QtConcurrent>
#include <QMessageBox>
#include <QTextEdit>
#include <QTableWidget>
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

public slots:
    void jobFinished(Data result);
    void printOutput(QString string, QColor c);
    void greyOutput();
    void enableRun();
    void disableRun();
    void stopJob();

private slots:

    //-----------
    //MainTab.cpp
    //-----------

    Request createRequest();

    void on_main_vcfDirBtn_clicked();
    void on_main_sampleDirBtn_clicked();
    void on_main_bedDirBtn_clicked();
    void on_main_runBtn_clicked();

    void on_main_testRareCastBtn_toggled(bool checked);
    void on_main_testRareSkatBtn_toggled(bool checked);
    void on_main_testBootChk_stateChanged(int arg1);
    void on_main_testBootChk_toggled(bool checked);
    void on_main_vcfWholeFileChk_toggled(bool checked);

    //-----------
    //SimulationTab.cpp
    //-----------
    void simulationTabInit();
    std::vector<SimulationRequestGroup> constructGroups(int nsteps, int highLow, Family family);

    void on_sim_runBtn_clicked();
    void on_sim_stopBtn_clicked();

    void addGroup(QTableWidget* table, QString n, QString cohort, QString depth, QString sdDepth, QString errorRate);
    void on_sim_groupAddBtn_clicked();
    void on_sim_groupRemoveBtn_clicked();

    void simEnableRare(bool valid);
    void on_sim_testRareCastBtn_toggled(bool checked);
    void on_sim_testRareSkatBtn_toggled(bool checked);

    void simulationFinished(Data results, SimulationRequest reqs);

    void on_main_stopBtn_clicked();

    void on_main_bedCollapseKBtn_toggled(bool checked);
    void on_main_bedCollapseExonBtn_toggled(bool checked);
    void on_main_bedCollapseGeneBtn_toggled(bool checked);

    void on_pushButton_pressed();

    void switchToCaseControl();
    void switchToQuantitative();
    void on_sim_caseControlBtn_pressed();
    void on_sim_quantitativeBtn_pressed();
    void on_sim_simulationSld_valueChanged(int value);

signals:
    void sendPlotData(QVector<double> values);

private:
    Ui::MainWindow *ui;
    QThread* jobThread;

    int getNumberOfStepsSim();
    int getHighLowCutoffSim();
    Family getFamilySim();
    int getNumberOfThreadsSim();
    Statistic getTestSim();
    int getBootstrapIterationsSim();
    int getCollapseSizeSim();
    bool getEarlyStopSim();

    QColor green = QColor::fromRgb(82, 145, 87);
    QColor grey = QColor::fromRgb(204, 205, 209);
    int plotCount = 1;

    QString prevOddsRatio = "1.0";
    QString prevR2 = "0.0";
    QVector<QString> prevCaseControlValues;
};

#endif // MAINWINDOW_H
