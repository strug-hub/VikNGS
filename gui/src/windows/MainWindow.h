#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "../src/vikNGS.h"
#include "../src/Variant.h"
#include "../simulation/simulation.h"
#include "../AsyncJob.h"

#include <QMainWindow>
#include <QFileDialog>
#include <QThread>
#include <QtConcurrent/QtConcurrent>
#include <QMessageBox>
#include <QTextEdit>
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
    void on_main_testRareCalphaBtn_toggled(bool checked);
    void on_main_testBootChk_stateChanged(int arg1);
    void on_main_testBootChk_toggled(bool checked);
    void on_main_vcfWholeFileChk_toggled(bool checked);

    void on_main_randomBtn_pressed();

    //-----------
    //SimulationTab.cpp
    //-----------
    void simulationTabInit();
    std::vector<SimulationRequestGroup> constructGroups(int ntests);
    SimulationRequest constructRequest(std::vector<SimulationRequestGroup> groups);

    void on_sim_runBtn_clicked();
    void on_sim_stopBtn_clicked();

    void addGroup(QString n, bool control, QString depth, QString sdDepth, QString errorRate);
    void on_sim_groupAddBtn_clicked();
    void on_sim_groupRemoveBtn_clicked();

    void simEnableRare(bool valid);
    void on_sim_testRareCastBtn_toggled(bool checked);
    void on_sim_testRareCalphaBtn_toggled(bool checked);

    void simulationFinished(std::vector<std::vector<Variant>> variants, SimulationRequest reqs);

    //-----------
    //qSimulationTab.cpp
    //-----------

    std::vector<SimulationRequestGroup> qConstructGroups(int run, int ntests);
    SimulationRequest qConstructRequest(std::vector<SimulationRequestGroup> groups);
    void qAddGroup(QString n, bool control, QString depth, QString sdDepth, QString errorRate);

    void qSimEnableRare(bool value);

    void on_qsim_runBtn_clicked();
    void on_qsim_groupAddBtn_clicked();
    void on_qsim_groupRemoveBtn_clicked();

    void on_qsim_stopBtn_clicked();

    void on_qsim_testRareCastBtn_toggled(bool checked);
    void on_qsim_testRareCalphaBtn_toggled(bool checked);

signals:
    void sendPlotData(QVector<double> values);

private:
    Ui::MainWindow *ui;
    QThread* jobThread;

    QColor green = QColor::fromRgb(82, 145, 87);
    QColor grey = QColor::fromRgb(204, 205, 209);
    int plotCount = 1;

    void warningDialog(QString message){

        QMessageBox *dialog = new QMessageBox;
        dialog->setWindowTitle("Warning");
        dialog->setIcon(QMessageBox::Warning);
        dialog->setText(message);
        dialog->show();
        return;

    }
};

#endif // MAINWINDOW_H
