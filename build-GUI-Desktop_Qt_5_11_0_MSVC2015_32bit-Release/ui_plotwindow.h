/********************************************************************************
** Form generated from reading UI file 'plotwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.11.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PLOTWINDOW_H
#define UI_PLOTWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "src/widgets/qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_PlotWindow
{
public:
    QVBoxLayout *verticalLayout;
    QCustomPlot *plot_chrPlt;
    QCustomPlot *plot_genomePlt;
    QHBoxLayout *GroupLayout;
    QGroupBox *plot_infoGrp;
    QVBoxLayout *verticalLayout_3;
    QLabel *plot_variantInfoLbl;
    QHBoxLayout *horizontalLayout_5;
    QHBoxLayout *horizontalLayout_2;
    QLabel *plot_altLbl;
    QLineEdit *plot_altTxt;
    QHBoxLayout *horizontalLayout_4;
    QLabel *plot_refLbl;
    QLineEdit *plot_refTxt;
    QGroupBox *plot_secondGrp;
    QGridLayout *gridLayout_2;
    QLabel *plot_title;

    void setupUi(QWidget *PlotWindow)
    {
        if (PlotWindow->objectName().isEmpty())
            PlotWindow->setObjectName(QStringLiteral("PlotWindow"));
        PlotWindow->resize(1171, 762);
        verticalLayout = new QVBoxLayout(PlotWindow);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        plot_chrPlt = new QCustomPlot(PlotWindow);
        plot_chrPlt->setObjectName(QStringLiteral("plot_chrPlt"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(85);
        sizePolicy.setHeightForWidth(plot_chrPlt->sizePolicy().hasHeightForWidth());
        plot_chrPlt->setSizePolicy(sizePolicy);

        verticalLayout->addWidget(plot_chrPlt);

        plot_genomePlt = new QCustomPlot(PlotWindow);
        plot_genomePlt->setObjectName(QStringLiteral("plot_genomePlt"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(10);
        sizePolicy1.setHeightForWidth(plot_genomePlt->sizePolicy().hasHeightForWidth());
        plot_genomePlt->setSizePolicy(sizePolicy1);

        verticalLayout->addWidget(plot_genomePlt);

        GroupLayout = new QHBoxLayout();
        GroupLayout->setObjectName(QStringLiteral("GroupLayout"));
        plot_infoGrp = new QGroupBox(PlotWindow);
        plot_infoGrp->setObjectName(QStringLiteral("plot_infoGrp"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(plot_infoGrp->sizePolicy().hasHeightForWidth());
        plot_infoGrp->setSizePolicy(sizePolicy2);
        verticalLayout_3 = new QVBoxLayout(plot_infoGrp);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        plot_variantInfoLbl = new QLabel(plot_infoGrp);
        plot_variantInfoLbl->setObjectName(QStringLiteral("plot_variantInfoLbl"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(1);
        sizePolicy3.setHeightForWidth(plot_variantInfoLbl->sizePolicy().hasHeightForWidth());
        plot_variantInfoLbl->setSizePolicy(sizePolicy3);
        QFont font;
        font.setPointSize(25);
        plot_variantInfoLbl->setFont(font);

        verticalLayout_3->addWidget(plot_variantInfoLbl);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        plot_altLbl = new QLabel(plot_infoGrp);
        plot_altLbl->setObjectName(QStringLiteral("plot_altLbl"));

        horizontalLayout_2->addWidget(plot_altLbl);

        plot_altTxt = new QLineEdit(plot_infoGrp);
        plot_altTxt->setObjectName(QStringLiteral("plot_altTxt"));
        plot_altTxt->setReadOnly(true);

        horizontalLayout_2->addWidget(plot_altTxt);


        horizontalLayout_5->addLayout(horizontalLayout_2);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        plot_refLbl = new QLabel(plot_infoGrp);
        plot_refLbl->setObjectName(QStringLiteral("plot_refLbl"));

        horizontalLayout_4->addWidget(plot_refLbl);

        plot_refTxt = new QLineEdit(plot_infoGrp);
        plot_refTxt->setObjectName(QStringLiteral("plot_refTxt"));

        horizontalLayout_4->addWidget(plot_refTxt);


        horizontalLayout_5->addLayout(horizontalLayout_4);


        verticalLayout_3->addLayout(horizontalLayout_5);


        GroupLayout->addWidget(plot_infoGrp);

        plot_secondGrp = new QGroupBox(PlotWindow);
        plot_secondGrp->setObjectName(QStringLiteral("plot_secondGrp"));
        sizePolicy2.setHeightForWidth(plot_secondGrp->sizePolicy().hasHeightForWidth());
        plot_secondGrp->setSizePolicy(sizePolicy2);
        gridLayout_2 = new QGridLayout(plot_secondGrp);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        plot_title = new QLabel(plot_secondGrp);
        plot_title->setObjectName(QStringLiteral("plot_title"));
        sizePolicy3.setHeightForWidth(plot_title->sizePolicy().hasHeightForWidth());
        plot_title->setSizePolicy(sizePolicy3);
        plot_title->setFont(font);

        gridLayout_2->addWidget(plot_title, 0, 0, 1, 1);


        GroupLayout->addWidget(plot_secondGrp);

        GroupLayout->setStretch(0, 1);
        GroupLayout->setStretch(1, 1);

        verticalLayout->addLayout(GroupLayout);

        verticalLayout->setStretch(0, 65);
        verticalLayout->setStretch(1, 25);
        verticalLayout->setStretch(2, 8);

        retranslateUi(PlotWindow);

        QMetaObject::connectSlotsByName(PlotWindow);
    } // setupUi

    void retranslateUi(QWidget *PlotWindow)
    {
        PlotWindow->setWindowTitle(QApplication::translate("PlotWindow", "Plotter", nullptr));
        plot_infoGrp->setTitle(QApplication::translate("PlotWindow", "Variant Information", nullptr));
        plot_variantInfoLbl->setText(QString());
        plot_altLbl->setText(QApplication::translate("PlotWindow", "ALT", nullptr));
        plot_refLbl->setText(QApplication::translate("PlotWindow", "REF", nullptr));
        plot_secondGrp->setTitle(QApplication::translate("PlotWindow", "Group", nullptr));
        plot_title->setText(QApplication::translate("PlotWindow", "Title", nullptr));
    } // retranslateUi

};

namespace Ui {
    class PlotWindow: public Ui_PlotWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PLOTWINDOW_H
