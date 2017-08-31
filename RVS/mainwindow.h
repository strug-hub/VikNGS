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

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
