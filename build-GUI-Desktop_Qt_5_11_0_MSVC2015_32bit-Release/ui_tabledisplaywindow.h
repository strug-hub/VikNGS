/********************************************************************************
** Form generated from reading UI file 'tabledisplaywindow.ui'
**
** Created by: Qt User Interface Compiler version 5.11.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_TABLEDISPLAYWINDOW_H
#define UI_TABLEDISPLAYWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_TableDisplayWindow
{
public:
    QHBoxLayout *horizontalLayout;
    QGroupBox *table_variantGrp;
    QVBoxLayout *verticalLayout_5;
    QTableWidget *table_variantTbl;
    QGroupBox *table_genoGrp;
    QVBoxLayout *verticalLayout;
    QTableWidget *table_genoTbl;
    QGroupBox *table_genoFreqGrp;
    QVBoxLayout *verticalLayout_4;
    QTableWidget *table_genoFreqTbl;

    void setupUi(QWidget *TableDisplayWindow)
    {
        if (TableDisplayWindow->objectName().isEmpty())
            TableDisplayWindow->setObjectName(QStringLiteral("TableDisplayWindow"));
        TableDisplayWindow->resize(721, 575);
        horizontalLayout = new QHBoxLayout(TableDisplayWindow);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        table_variantGrp = new QGroupBox(TableDisplayWindow);
        table_variantGrp->setObjectName(QStringLiteral("table_variantGrp"));
        verticalLayout_5 = new QVBoxLayout(table_variantGrp);
        verticalLayout_5->setObjectName(QStringLiteral("verticalLayout_5"));
        table_variantTbl = new QTableWidget(table_variantGrp);
        if (table_variantTbl->columnCount() < 4)
            table_variantTbl->setColumnCount(4);
        QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
        table_variantTbl->setHorizontalHeaderItem(0, __qtablewidgetitem);
        QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
        table_variantTbl->setHorizontalHeaderItem(1, __qtablewidgetitem1);
        QTableWidgetItem *__qtablewidgetitem2 = new QTableWidgetItem();
        table_variantTbl->setHorizontalHeaderItem(2, __qtablewidgetitem2);
        QTableWidgetItem *__qtablewidgetitem3 = new QTableWidgetItem();
        table_variantTbl->setHorizontalHeaderItem(3, __qtablewidgetitem3);
        table_variantTbl->setObjectName(QStringLiteral("table_variantTbl"));
        table_variantTbl->setSelectionBehavior(QAbstractItemView::SelectRows);

        verticalLayout_5->addWidget(table_variantTbl);


        horizontalLayout->addWidget(table_variantGrp);

        table_genoGrp = new QGroupBox(TableDisplayWindow);
        table_genoGrp->setObjectName(QStringLiteral("table_genoGrp"));
        verticalLayout = new QVBoxLayout(table_genoGrp);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        table_genoTbl = new QTableWidget(table_genoGrp);
        if (table_genoTbl->columnCount() < 4)
            table_genoTbl->setColumnCount(4);
        QTableWidgetItem *__qtablewidgetitem4 = new QTableWidgetItem();
        table_genoTbl->setHorizontalHeaderItem(0, __qtablewidgetitem4);
        QTableWidgetItem *__qtablewidgetitem5 = new QTableWidgetItem();
        table_genoTbl->setHorizontalHeaderItem(1, __qtablewidgetitem5);
        QTableWidgetItem *__qtablewidgetitem6 = new QTableWidgetItem();
        table_genoTbl->setHorizontalHeaderItem(2, __qtablewidgetitem6);
        QTableWidgetItem *__qtablewidgetitem7 = new QTableWidgetItem();
        table_genoTbl->setHorizontalHeaderItem(3, __qtablewidgetitem7);
        table_genoTbl->setObjectName(QStringLiteral("table_genoTbl"));
        table_genoTbl->setSelectionBehavior(QAbstractItemView::SelectRows);

        verticalLayout->addWidget(table_genoTbl);


        horizontalLayout->addWidget(table_genoGrp);

        table_genoFreqGrp = new QGroupBox(TableDisplayWindow);
        table_genoFreqGrp->setObjectName(QStringLiteral("table_genoFreqGrp"));
        verticalLayout_4 = new QVBoxLayout(table_genoFreqGrp);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        table_genoFreqTbl = new QTableWidget(table_genoFreqGrp);
        if (table_genoFreqTbl->columnCount() < 4)
            table_genoFreqTbl->setColumnCount(4);
        QTableWidgetItem *__qtablewidgetitem8 = new QTableWidgetItem();
        table_genoFreqTbl->setHorizontalHeaderItem(0, __qtablewidgetitem8);
        QTableWidgetItem *__qtablewidgetitem9 = new QTableWidgetItem();
        table_genoFreqTbl->setHorizontalHeaderItem(1, __qtablewidgetitem9);
        QTableWidgetItem *__qtablewidgetitem10 = new QTableWidgetItem();
        table_genoFreqTbl->setHorizontalHeaderItem(2, __qtablewidgetitem10);
        QTableWidgetItem *__qtablewidgetitem11 = new QTableWidgetItem();
        table_genoFreqTbl->setHorizontalHeaderItem(3, __qtablewidgetitem11);
        table_genoFreqTbl->setObjectName(QStringLiteral("table_genoFreqTbl"));

        verticalLayout_4->addWidget(table_genoFreqTbl);


        horizontalLayout->addWidget(table_genoFreqGrp);


        retranslateUi(TableDisplayWindow);

        QMetaObject::connectSlotsByName(TableDisplayWindow);
    } // setupUi

    void retranslateUi(QWidget *TableDisplayWindow)
    {
        TableDisplayWindow->setWindowTitle(QApplication::translate("TableDisplayWindow", "Form", nullptr));
        table_variantGrp->setTitle(QApplication::translate("TableDisplayWindow", "Variant - Genotype Frequency", nullptr));
        QTableWidgetItem *___qtablewidgetitem = table_variantTbl->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QApplication::translate("TableDisplayWindow", "MAF", nullptr));
        QTableWidgetItem *___qtablewidgetitem1 = table_variantTbl->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QApplication::translate("TableDisplayWindow", "0", nullptr));
        QTableWidgetItem *___qtablewidgetitem2 = table_variantTbl->horizontalHeaderItem(2);
        ___qtablewidgetitem2->setText(QApplication::translate("TableDisplayWindow", "1", nullptr));
        QTableWidgetItem *___qtablewidgetitem3 = table_variantTbl->horizontalHeaderItem(3);
        ___qtablewidgetitem3->setText(QApplication::translate("TableDisplayWindow", "2", nullptr));
        table_genoGrp->setTitle(QApplication::translate("TableDisplayWindow", "Genotypes", nullptr));
        QTableWidgetItem *___qtablewidgetitem4 = table_genoTbl->horizontalHeaderItem(0);
        ___qtablewidgetitem4->setText(QApplication::translate("TableDisplayWindow", "Sample", nullptr));
        QTableWidgetItem *___qtablewidgetitem5 = table_genoTbl->horizontalHeaderItem(1);
        ___qtablewidgetitem5->setText(QApplication::translate("TableDisplayWindow", "True GT", nullptr));
        QTableWidgetItem *___qtablewidgetitem6 = table_genoTbl->horizontalHeaderItem(2);
        ___qtablewidgetitem6->setText(QApplication::translate("TableDisplayWindow", "GT Call", nullptr));
        QTableWidgetItem *___qtablewidgetitem7 = table_genoTbl->horizontalHeaderItem(3);
        ___qtablewidgetitem7->setText(QApplication::translate("TableDisplayWindow", "Expected GT", nullptr));
        table_genoFreqGrp->setTitle(QApplication::translate("TableDisplayWindow", "Genotype Frequency", nullptr));
        QTableWidgetItem *___qtablewidgetitem8 = table_genoFreqTbl->horizontalHeaderItem(0);
        ___qtablewidgetitem8->setText(QApplication::translate("TableDisplayWindow", "0", nullptr));
        QTableWidgetItem *___qtablewidgetitem9 = table_genoFreqTbl->horizontalHeaderItem(1);
        ___qtablewidgetitem9->setText(QApplication::translate("TableDisplayWindow", "1", nullptr));
        QTableWidgetItem *___qtablewidgetitem10 = table_genoFreqTbl->horizontalHeaderItem(2);
        ___qtablewidgetitem10->setText(QApplication::translate("TableDisplayWindow", "2", nullptr));
        QTableWidgetItem *___qtablewidgetitem11 = table_genoFreqTbl->horizontalHeaderItem(3);
        ___qtablewidgetitem11->setText(QApplication::translate("TableDisplayWindow", "Average", nullptr));
    } // retranslateUi

};

namespace Ui {
    class TableDisplayWindow: public Ui_TableDisplayWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TABLEDISPLAYWINDOW_H
