#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "RVS.h"

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
    void on_vcfBrowse_clicked();

    void on_runSimulationButton_clicked();

    void on_addGroupButton_clicked();

    void on_removeGroupButton_clicked();

    void on_binSlider_valueChanged(int value);

    void replot();

    void on_rareRdoButton_clicked();

    void on_bootstrapChkBox_stateChanged(int arg1);

    void on_sampleBrowse_clicked();

    void on_sampleBrowse_2_clicked();

    void on_runButton_clicked();

    void replotMain();

    void on_mainBinSlider_valueChanged(int value);


private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
