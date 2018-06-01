#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "../src/RVS.h"
#include "../src/Variant.h"
#include "../simulation/simulation.h"
#include "../runner.h"

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
    void jobFinished(QVector<Variant> variants);
    void printOutput(QString string, QColor c);
    void greyOutput();


private slots:
    double calculatePower(QVector<double> pval);
    std::vector<SimulationRequestGroup> constructGroups(int run, int ntests);
    SimulationRequest constructRequest(std::vector<SimulationRequestGroup> groups);

    void on_sim_runBtn_clicked();
    void on_sim_groupAddBtn_clicked();
    void on_sim_groupRemoveBtn_clicked();

    void on_sim_testBootChk_stateChanged(int arg1);
    void on_sim_testBootChk_toggled(bool checked);
    void on_sim_testRareCastBtn_toggled(bool checked);
    void on_sim_testRareCalphaBtn_toggled(bool checked);

    void simulationFinished(std::vector<std::vector<Variant>> variants, std::vector<SimulationRequest> reqs);

    //--------------
    void on_main_vcfDirBtn_clicked();
    void on_main_sampleDirBtn_clicked();
    void on_main_bedDirBtn_clicked();

    void on_main_runBtn_clicked();

    void on_main_testRareCastBtn_toggled(bool checked);
    void on_main_testRareCalphaBtn_toggled(bool checked);
    void on_main_testBootChk_stateChanged(int arg1);
    void on_main_testBootChk_toggled(bool checked);
    void on_main_vcfWholeFileChk_toggled(bool checked);

signals:
    void sendPlotData(QVector<double> values);

private:
    Ui::MainWindow *ui;
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
