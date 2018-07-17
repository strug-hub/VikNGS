/********************************************************************************
** Form generated from reading UI file 'simplotwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.11.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SIMPLOTWINDOW_H
#define UI_SIMPLOTWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDial>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLCDNumber>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "src/widgets/qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_SimPlotWindow
{
public:
    QHBoxLayout *horizontalLayout_8;
    QVBoxLayout *verticalLayout_3;
    QHBoxLayout *horizontalLayout;
    QCustomPlot *simplot_plot1;
    QCustomPlot *simplot_plot2;
    QGroupBox *simplot_legendGrp;
    QGridLayout *gridLayout_2;
    QLabel *simplot_legend1Lbl;
    QLabel *simplot_legend2Lbl;
    QWidget *simplot_legend4Sqr;
    QWidget *simplot_legend1Sqr;
    QWidget *simplot_legend3Sqr;
    QWidget *simplot_legend2Sqr;
    QLabel *simplot_legend3Lbl;
    QLabel *simplot_legend4Lbl;
    QHBoxLayout *horizontalLayout_4;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout_2;
    QDial *simplot_alphaDial;
    QHBoxLayout *horizontalLayout_7;
    QLineEdit *simplot_alphaTxt;
    QGroupBox *simplot_measureGrp;
    QVBoxLayout *verticalLayout_4;
    QTableWidget *simplot_powerTbl;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout;
    QLabel *simplot_ncasesLbl;
    QLabel *simplot_ncontrolsLbl;
    QLCDNumber *simplot_ncontrolsDgt;
    QLCDNumber *simplot_ncasesDgt;
    QGroupBox *simplot_genoGrp;
    QVBoxLayout *verticalLayout;
    QPushButton *pushButton;

    void setupUi(QWidget *SimPlotWindow)
    {
        if (SimPlotWindow->objectName().isEmpty())
            SimPlotWindow->setObjectName(QStringLiteral("SimPlotWindow"));
        SimPlotWindow->resize(1102, 691);
        horizontalLayout_8 = new QHBoxLayout(SimPlotWindow);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        simplot_plot1 = new QCustomPlot(SimPlotWindow);
        simplot_plot1->setObjectName(QStringLiteral("simplot_plot1"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(85);
        sizePolicy.setHeightForWidth(simplot_plot1->sizePolicy().hasHeightForWidth());
        simplot_plot1->setSizePolicy(sizePolicy);

        horizontalLayout->addWidget(simplot_plot1);

        simplot_plot2 = new QCustomPlot(SimPlotWindow);
        simplot_plot2->setObjectName(QStringLiteral("simplot_plot2"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(10);
        sizePolicy1.setHeightForWidth(simplot_plot2->sizePolicy().hasHeightForWidth());
        simplot_plot2->setSizePolicy(sizePolicy1);

        horizontalLayout->addWidget(simplot_plot2);


        verticalLayout_3->addLayout(horizontalLayout);

        simplot_legendGrp = new QGroupBox(SimPlotWindow);
        simplot_legendGrp->setObjectName(QStringLiteral("simplot_legendGrp"));
        gridLayout_2 = new QGridLayout(simplot_legendGrp);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        simplot_legend1Lbl = new QLabel(simplot_legendGrp);
        simplot_legend1Lbl->setObjectName(QStringLiteral("simplot_legend1Lbl"));
        simplot_legend1Lbl->setMinimumSize(QSize(24, 24));

        gridLayout_2->addWidget(simplot_legend1Lbl, 0, 1, 1, 1);

        simplot_legend2Lbl = new QLabel(simplot_legendGrp);
        simplot_legend2Lbl->setObjectName(QStringLiteral("simplot_legend2Lbl"));

        gridLayout_2->addWidget(simplot_legend2Lbl, 0, 3, 1, 1);

        simplot_legend4Sqr = new QWidget(simplot_legendGrp);
        simplot_legend4Sqr->setObjectName(QStringLiteral("simplot_legend4Sqr"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(1);
        sizePolicy2.setVerticalStretch(1);
        sizePolicy2.setHeightForWidth(simplot_legend4Sqr->sizePolicy().hasHeightForWidth());
        simplot_legend4Sqr->setSizePolicy(sizePolicy2);
        simplot_legend4Sqr->setMinimumSize(QSize(24, 24));
        simplot_legend4Sqr->setMaximumSize(QSize(48, 48));
        simplot_legend4Sqr->setSizeIncrement(QSize(0, 0));

        gridLayout_2->addWidget(simplot_legend4Sqr, 0, 6, 1, 1);

        simplot_legend1Sqr = new QWidget(simplot_legendGrp);
        simplot_legend1Sqr->setObjectName(QStringLiteral("simplot_legend1Sqr"));
        sizePolicy2.setHeightForWidth(simplot_legend1Sqr->sizePolicy().hasHeightForWidth());
        simplot_legend1Sqr->setSizePolicy(sizePolicy2);
        simplot_legend1Sqr->setMinimumSize(QSize(24, 24));
        simplot_legend1Sqr->setMaximumSize(QSize(48, 48));
        simplot_legend1Sqr->setSizeIncrement(QSize(0, 0));

        gridLayout_2->addWidget(simplot_legend1Sqr, 0, 0, 1, 1);

        simplot_legend3Sqr = new QWidget(simplot_legendGrp);
        simplot_legend3Sqr->setObjectName(QStringLiteral("simplot_legend3Sqr"));
        sizePolicy2.setHeightForWidth(simplot_legend3Sqr->sizePolicy().hasHeightForWidth());
        simplot_legend3Sqr->setSizePolicy(sizePolicy2);
        simplot_legend3Sqr->setMinimumSize(QSize(24, 24));
        simplot_legend3Sqr->setMaximumSize(QSize(48, 48));
        simplot_legend3Sqr->setSizeIncrement(QSize(0, 0));

        gridLayout_2->addWidget(simplot_legend3Sqr, 0, 4, 1, 1);

        simplot_legend2Sqr = new QWidget(simplot_legendGrp);
        simplot_legend2Sqr->setObjectName(QStringLiteral("simplot_legend2Sqr"));
        sizePolicy2.setHeightForWidth(simplot_legend2Sqr->sizePolicy().hasHeightForWidth());
        simplot_legend2Sqr->setSizePolicy(sizePolicy2);
        simplot_legend2Sqr->setMinimumSize(QSize(24, 24));
        simplot_legend2Sqr->setMaximumSize(QSize(48, 48));
        simplot_legend2Sqr->setSizeIncrement(QSize(0, 0));

        gridLayout_2->addWidget(simplot_legend2Sqr, 0, 2, 1, 1);

        simplot_legend3Lbl = new QLabel(simplot_legendGrp);
        simplot_legend3Lbl->setObjectName(QStringLiteral("simplot_legend3Lbl"));

        gridLayout_2->addWidget(simplot_legend3Lbl, 0, 5, 1, 1);

        simplot_legend4Lbl = new QLabel(simplot_legendGrp);
        simplot_legend4Lbl->setObjectName(QStringLiteral("simplot_legend4Lbl"));

        gridLayout_2->addWidget(simplot_legend4Lbl, 0, 7, 1, 1);


        verticalLayout_3->addWidget(simplot_legendGrp);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        groupBox = new QGroupBox(SimPlotWindow);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        verticalLayout_2 = new QVBoxLayout(groupBox);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        simplot_alphaDial = new QDial(groupBox);
        simplot_alphaDial->setObjectName(QStringLiteral("simplot_alphaDial"));
        simplot_alphaDial->setMinimum(1);
        simplot_alphaDial->setMaximum(10000);
        simplot_alphaDial->setValue(1301);
        simplot_alphaDial->setTracking(true);
        simplot_alphaDial->setOrientation(Qt::Horizontal);
        simplot_alphaDial->setInvertedAppearance(false);
        simplot_alphaDial->setWrapping(false);
        simplot_alphaDial->setNotchTarget(100);
        simplot_alphaDial->setNotchesVisible(true);

        verticalLayout_2->addWidget(simplot_alphaDial);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        simplot_alphaTxt = new QLineEdit(groupBox);
        simplot_alphaTxt->setObjectName(QStringLiteral("simplot_alphaTxt"));

        horizontalLayout_7->addWidget(simplot_alphaTxt);


        verticalLayout_2->addLayout(horizontalLayout_7);


        horizontalLayout_4->addWidget(groupBox);

        simplot_measureGrp = new QGroupBox(SimPlotWindow);
        simplot_measureGrp->setObjectName(QStringLiteral("simplot_measureGrp"));
        verticalLayout_4 = new QVBoxLayout(simplot_measureGrp);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        simplot_powerTbl = new QTableWidget(simplot_measureGrp);
        simplot_powerTbl->setObjectName(QStringLiteral("simplot_powerTbl"));
        simplot_powerTbl->setFrameShape(QFrame::NoFrame);
        simplot_powerTbl->setFrameShadow(QFrame::Plain);
        simplot_powerTbl->setSelectionMode(QAbstractItemView::NoSelection);
        simplot_powerTbl->setSelectionBehavior(QAbstractItemView::SelectRows);
        simplot_powerTbl->setCornerButtonEnabled(true);
        simplot_powerTbl->horizontalHeader()->setVisible(false);
        simplot_powerTbl->verticalHeader()->setVisible(false);
        simplot_powerTbl->verticalHeader()->setStretchLastSection(false);

        verticalLayout_4->addWidget(simplot_powerTbl);


        horizontalLayout_4->addWidget(simplot_measureGrp);

        groupBox_2 = new QGroupBox(SimPlotWindow);
        groupBox_2->setObjectName(QStringLiteral("groupBox_2"));
        gridLayout = new QGridLayout(groupBox_2);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        simplot_ncasesLbl = new QLabel(groupBox_2);
        simplot_ncasesLbl->setObjectName(QStringLiteral("simplot_ncasesLbl"));
        simplot_ncasesLbl->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(simplot_ncasesLbl, 2, 0, 1, 1);

        simplot_ncontrolsLbl = new QLabel(groupBox_2);
        simplot_ncontrolsLbl->setObjectName(QStringLiteral("simplot_ncontrolsLbl"));
        simplot_ncontrolsLbl->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(simplot_ncontrolsLbl, 4, 0, 1, 1);

        simplot_ncontrolsDgt = new QLCDNumber(groupBox_2);
        simplot_ncontrolsDgt->setObjectName(QStringLiteral("simplot_ncontrolsDgt"));
        simplot_ncontrolsDgt->setFrameShape(QFrame::NoFrame);
        simplot_ncontrolsDgt->setFrameShadow(QFrame::Sunken);
        simplot_ncontrolsDgt->setSmallDecimalPoint(false);
        simplot_ncontrolsDgt->setSegmentStyle(QLCDNumber::Flat);

        gridLayout->addWidget(simplot_ncontrolsDgt, 3, 0, 1, 1);

        simplot_ncasesDgt = new QLCDNumber(groupBox_2);
        simplot_ncasesDgt->setObjectName(QStringLiteral("simplot_ncasesDgt"));
        simplot_ncasesDgt->setFrameShape(QFrame::NoFrame);
        simplot_ncasesDgt->setSegmentStyle(QLCDNumber::Flat);

        gridLayout->addWidget(simplot_ncasesDgt, 0, 0, 1, 1);


        horizontalLayout_4->addWidget(groupBox_2);

        simplot_genoGrp = new QGroupBox(SimPlotWindow);
        simplot_genoGrp->setObjectName(QStringLiteral("simplot_genoGrp"));
        verticalLayout = new QVBoxLayout(simplot_genoGrp);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        pushButton = new QPushButton(simplot_genoGrp);
        pushButton->setObjectName(QStringLiteral("pushButton"));

        verticalLayout->addWidget(pushButton);


        horizontalLayout_4->addWidget(simplot_genoGrp);

        horizontalLayout_4->setStretch(0, 1);
        horizontalLayout_4->setStretch(1, 1);
        horizontalLayout_4->setStretch(2, 1);
        horizontalLayout_4->setStretch(3, 3);

        verticalLayout_3->addLayout(horizontalLayout_4);

        verticalLayout_3->setStretch(0, 12);
        verticalLayout_3->setStretch(1, 1);
        verticalLayout_3->setStretch(2, 5);

        horizontalLayout_8->addLayout(verticalLayout_3);


        retranslateUi(SimPlotWindow);

        QMetaObject::connectSlotsByName(SimPlotWindow);
    } // setupUi

    void retranslateUi(QWidget *SimPlotWindow)
    {
        SimPlotWindow->setWindowTitle(QApplication::translate("SimPlotWindow", "Plotter", nullptr));
        simplot_legendGrp->setTitle(QApplication::translate("SimPlotWindow", "Legend", nullptr));
        simplot_legend1Lbl->setText(QApplication::translate("SimPlotWindow", "legend1", nullptr));
        simplot_legend2Lbl->setText(QApplication::translate("SimPlotWindow", "legend2", nullptr));
        simplot_legend3Lbl->setText(QApplication::translate("SimPlotWindow", "legend3", nullptr));
        simplot_legend4Lbl->setText(QApplication::translate("SimPlotWindow", "legend4", nullptr));
        groupBox->setTitle(QApplication::translate("SimPlotWindow", "Alpha level", nullptr));
        simplot_alphaTxt->setText(QApplication::translate("SimPlotWindow", "0.05", nullptr));
        simplot_measureGrp->setTitle(QApplication::translate("SimPlotWindow", "Label", nullptr));
        groupBox_2->setTitle(QApplication::translate("SimPlotWindow", "Sample Size", nullptr));
        simplot_ncasesLbl->setText(QApplication::translate("SimPlotWindow", "Cases", nullptr));
        simplot_ncontrolsLbl->setText(QApplication::translate("SimPlotWindow", "Controls", nullptr));
        simplot_genoGrp->setTitle(QApplication::translate("SimPlotWindow", "Genotypes", nullptr));
        pushButton->setText(QApplication::translate("SimPlotWindow", "Explore Genotypes", nullptr));
    } // retranslateUi

};

namespace Ui {
    class SimPlotWindow: public Ui_SimPlotWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SIMPLOTWINDOW_H
