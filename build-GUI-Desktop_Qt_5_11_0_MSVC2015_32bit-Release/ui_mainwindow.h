/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.11.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *MainWindowFrame;
    QVBoxLayout *verticalLayout_8;
    QTabWidget *tabWidget;
    QWidget *mainTab;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *InputLayout1;
    QGroupBox *main_vcfGrp;
    QVBoxLayout *verticalLayout_7;
    QHBoxLayout *VCFDirLayout;
    QLineEdit *main_vcfDirTxt;
    QPushButton *main_vcfDirBtn;
    QGridLayout *gridLayout;
    QLabel *main_vcfMissingLbl;
    QLineEdit *main_vcfFromPosTxt;
    QCheckBox *main_vcfWholeFileChk;
    QLineEdit *main_vcfChromFilterTxt;
    QLineEdit *main_vcfToPosTxt;
    QLabel *main_vcfToPosLbl;
    QLineEdit *main_vcfMissingTxt;
    QCheckBox *main_vcfPassChk;
    QLabel *main_vcfFromPosLbl;
    QLabel *main_vcfMafLbl;
    QLineEdit *main_vcfMafTxt;
    QLabel *main_vcfChrLbl;
    QGroupBox *main_sampleGrp;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *SampleDirLayout;
    QLineEdit *main_sampleDirTxt;
    QPushButton *main_sampleDirBtn;
    QHBoxLayout *SampleParamLayout;
    QLineEdit *main_sampleDepthTxt;
    QLabel *main_sampleDepthLbl;
    QSpacerItem *verticalSpacer_2;
    QSpacerItem *horizontalSpacer;
    QVBoxLayout *verticalLayout_5;
    QGroupBox *main_bedGrp;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *BEDDirLayout;
    QLineEdit *main_bedDirTxt;
    QPushButton *main_bedDirBtn;
    QGroupBox *main_bedCollapseGrp;
    QHBoxLayout *horizontalLayout_13;
    QHBoxLayout *horizontalLayout_5;
    QRadioButton *main_bedCollapseExonBtn;
    QRadioButton *main_bedCollapseGeneBtn;
    QRadioButton *main_bedCollapseKBtn;
    QLineEdit *main_bedCollapseKTxt;
    QGroupBox *main_testGrp;
    QHBoxLayout *horizontalLayout_12;
    QVBoxLayout *TestLayout;
    QRadioButton *main_testCommonBtn;
    QRadioButton *main_testRareCastBtn;
    QRadioButton *main_testRareCalphaBtn;
    QVBoxLayout *verticalLayout_4;
    QCheckBox *main_testBootChk;
    QHBoxLayout *horizontalLayout_7;
    QLineEdit *main_testBootTxt;
    QLabel *main_testBootLbl;
    QCheckBox *main_testStopChk;
    QGroupBox *main_startGrp;
    QHBoxLayout *horizontalLayout;
    QFormLayout *formLayout;
    QLineEdit *main_threadsTxt;
    QLabel *main_threadsLbl;
    QLineEdit *main_batchSizeTxt;
    QLabel *main_batchSizeLbl;
    QCheckBox *main_plotChk;
    QCheckBox *main_rvsChk;
    QCheckBox *main_gtChk;
    QPushButton *main_runBtn;
    QPushButton *main_stopBtn;
    QPushButton *main_randomBtn;
    QLineEdit *main_nrandomTxt;
    QWidget *simulationTab;
    QHBoxLayout *horizontalLayout_8;
    QGroupBox *sim_groupGrp;
    QVBoxLayout *verticalLayout_10;
    QVBoxLayout *verticalLayout_11;
    QPushButton *sim_groupAddBtn;
    QPushButton *sim_groupRemoveBtn;
    QTableWidget *sim_groupTbl;
    QVBoxLayout *verticalLayout_3;
    QGroupBox *sim_variantGrp;
    QGridLayout *gridLayout_2;
    QLineEdit *sim_variantMafMaxTxt;
    QLineEdit *sim_variantSizeTxt;
    QLabel *sim_variantMafMinLbl;
    QLabel *sim_variantSizeLbl;
    QLabel *sim_variantMafMaxLbl;
    QLineEdit *sim_variantMafMinTxt;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_6;
    QLabel *sim_powerStepLbl;
    QLabel *sim_groupHighLowLbl;
    QLineEdit *sim_powerStepTxt;
    QLineEdit *sim_groupHighLowTxt;
    QLabel *sim_oddsRatioLbl;
    QLineEdit *sim_oddsRatioTxt;
    QGroupBox *sim_testGrp;
    QHBoxLayout *horizontalLayout_4;
    QVBoxLayout *verticalLayout_6;
    QRadioButton *sim_testCommonBtn;
    QRadioButton *sim_testRareCastBtn;
    QRadioButton *sim_testRareCalphaBtn;
    QVBoxLayout *verticalLayout_9;
    QHBoxLayout *horizontalLayout_3;
    QGridLayout *gridLayout_4;
    QLabel *sim_testBootLbl;
    QLabel *sim_collapseLbl;
    QLineEdit *sim_testBootTxt;
    QLineEdit *sim_collapseTxt;
    QCheckBox *sim_testStopChk;
    QPushButton *sim_runBtn;
    QPushButton *sim_stopBtn;
    QWidget *qSimTab;
    QHBoxLayout *horizontalLayout_10;
    QGroupBox *qsim_groupGrp;
    QVBoxLayout *verticalLayout_12;
    QVBoxLayout *verticalLayout_13;
    QPushButton *qsim_groupAddBtn;
    QPushButton *qsim_groupRemoveBtn;
    QTableWidget *qsim_groupTbl;
    QVBoxLayout *verticalLayout_14;
    QGroupBox *qsim_variantGrp;
    QGridLayout *gridLayout_3;
    QLineEdit *qsim_variantMafMaxTxt;
    QLineEdit *qsim_variantSizeTxt;
    QLabel *qsim_variantMafMinLbl;
    QLabel *qsim_variantSizeLbl;
    QLabel *qsim_variantMafMaxLbl;
    QLineEdit *qsim_variantMafMinTxt;
    QGroupBox *qsim_runGrp;
    QGridLayout *gridLayout_7;
    QLabel *qsim_powerStepLbl;
    QLabel *qsim_groupHighLowLbl;
    QLineEdit *qsim_powerStepTxt;
    QLineEdit *qsim_groupHighLowTxt;
    QLabel *qsim_R2Lbl;
    QLineEdit *qsim_R2Txt;
    QGroupBox *qsim_testGrp;
    QHBoxLayout *horizontalLayout_6;
    QVBoxLayout *verticalLayout_15;
    QRadioButton *qsim_testCommonBtn;
    QRadioButton *qsim_testRareCastBtn;
    QRadioButton *qsim_testRareCalphaBtn;
    QVBoxLayout *verticalLayout_16;
    QHBoxLayout *horizontalLayout_9;
    QGridLayout *gridLayout_5;
    QLabel *qsim_testBootLbl;
    QLabel *qsim_collapseLbl;
    QLineEdit *qsim_testBootTxt;
    QLineEdit *qsim_collapseTxt;
    QCheckBox *qsim_testStopChk;
    QPushButton *qsim_runBtn;
    QPushButton *qsim_stopBtn;
    QTextEdit *outputBox;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(900, 752);
        MainWindowFrame = new QWidget(MainWindow);
        MainWindowFrame->setObjectName(QStringLiteral("MainWindowFrame"));
        verticalLayout_8 = new QVBoxLayout(MainWindowFrame);
        verticalLayout_8->setSpacing(6);
        verticalLayout_8->setContentsMargins(11, 11, 11, 11);
        verticalLayout_8->setObjectName(QStringLiteral("verticalLayout_8"));
        tabWidget = new QTabWidget(MainWindowFrame);
        tabWidget->setObjectName(QStringLiteral("tabWidget"));
        tabWidget->setEnabled(true);
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(60);
        sizePolicy.setHeightForWidth(tabWidget->sizePolicy().hasHeightForWidth());
        tabWidget->setSizePolicy(sizePolicy);
        tabWidget->setMovable(true);
        tabWidget->setTabBarAutoHide(false);
        mainTab = new QWidget();
        mainTab->setObjectName(QStringLiteral("mainTab"));
        horizontalLayout_2 = new QHBoxLayout(mainTab);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        InputLayout1 = new QVBoxLayout();
        InputLayout1->setSpacing(6);
        InputLayout1->setObjectName(QStringLiteral("InputLayout1"));
        main_vcfGrp = new QGroupBox(mainTab);
        main_vcfGrp->setObjectName(QStringLiteral("main_vcfGrp"));
        verticalLayout_7 = new QVBoxLayout(main_vcfGrp);
        verticalLayout_7->setSpacing(6);
        verticalLayout_7->setContentsMargins(11, 11, 11, 11);
        verticalLayout_7->setObjectName(QStringLiteral("verticalLayout_7"));
        VCFDirLayout = new QHBoxLayout();
        VCFDirLayout->setSpacing(6);
        VCFDirLayout->setObjectName(QStringLiteral("VCFDirLayout"));
        main_vcfDirTxt = new QLineEdit(main_vcfGrp);
        main_vcfDirTxt->setObjectName(QStringLiteral("main_vcfDirTxt"));
        main_vcfDirTxt->setReadOnly(true);

        VCFDirLayout->addWidget(main_vcfDirTxt);

        main_vcfDirBtn = new QPushButton(main_vcfGrp);
        main_vcfDirBtn->setObjectName(QStringLiteral("main_vcfDirBtn"));
        main_vcfDirBtn->setMaximumSize(QSize(20, 16777215));

        VCFDirLayout->addWidget(main_vcfDirBtn);


        verticalLayout_7->addLayout(VCFDirLayout);

        gridLayout = new QGridLayout();
        gridLayout->setSpacing(6);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        main_vcfMissingLbl = new QLabel(main_vcfGrp);
        main_vcfMissingLbl->setObjectName(QStringLiteral("main_vcfMissingLbl"));

        gridLayout->addWidget(main_vcfMissingLbl, 1, 1, 1, 1);

        main_vcfFromPosTxt = new QLineEdit(main_vcfGrp);
        main_vcfFromPosTxt->setObjectName(QStringLiteral("main_vcfFromPosTxt"));
        main_vcfFromPosTxt->setEnabled(false);

        gridLayout->addWidget(main_vcfFromPosTxt, 5, 0, 1, 1);

        main_vcfWholeFileChk = new QCheckBox(main_vcfGrp);
        main_vcfWholeFileChk->setObjectName(QStringLiteral("main_vcfWholeFileChk"));
        main_vcfWholeFileChk->setChecked(true);

        gridLayout->addWidget(main_vcfWholeFileChk, 3, 0, 1, 1);

        main_vcfChromFilterTxt = new QLineEdit(main_vcfGrp);
        main_vcfChromFilterTxt->setObjectName(QStringLiteral("main_vcfChromFilterTxt"));
        main_vcfChromFilterTxt->setEnabled(false);

        gridLayout->addWidget(main_vcfChromFilterTxt, 4, 0, 1, 1);

        main_vcfToPosTxt = new QLineEdit(main_vcfGrp);
        main_vcfToPosTxt->setObjectName(QStringLiteral("main_vcfToPosTxt"));
        main_vcfToPosTxt->setEnabled(false);

        gridLayout->addWidget(main_vcfToPosTxt, 6, 0, 1, 1);

        main_vcfToPosLbl = new QLabel(main_vcfGrp);
        main_vcfToPosLbl->setObjectName(QStringLiteral("main_vcfToPosLbl"));
        main_vcfToPosLbl->setEnabled(false);

        gridLayout->addWidget(main_vcfToPosLbl, 6, 1, 1, 1);

        main_vcfMissingTxt = new QLineEdit(main_vcfGrp);
        main_vcfMissingTxt->setObjectName(QStringLiteral("main_vcfMissingTxt"));

        gridLayout->addWidget(main_vcfMissingTxt, 1, 0, 1, 1);

        main_vcfPassChk = new QCheckBox(main_vcfGrp);
        main_vcfPassChk->setObjectName(QStringLiteral("main_vcfPassChk"));
        main_vcfPassChk->setChecked(true);

        gridLayout->addWidget(main_vcfPassChk, 7, 0, 1, 1);

        main_vcfFromPosLbl = new QLabel(main_vcfGrp);
        main_vcfFromPosLbl->setObjectName(QStringLiteral("main_vcfFromPosLbl"));
        main_vcfFromPosLbl->setEnabled(false);

        gridLayout->addWidget(main_vcfFromPosLbl, 5, 1, 1, 1);

        main_vcfMafLbl = new QLabel(main_vcfGrp);
        main_vcfMafLbl->setObjectName(QStringLiteral("main_vcfMafLbl"));

        gridLayout->addWidget(main_vcfMafLbl, 0, 1, 1, 1);

        main_vcfMafTxt = new QLineEdit(main_vcfGrp);
        main_vcfMafTxt->setObjectName(QStringLiteral("main_vcfMafTxt"));

        gridLayout->addWidget(main_vcfMafTxt, 0, 0, 1, 1);

        main_vcfChrLbl = new QLabel(main_vcfGrp);
        main_vcfChrLbl->setObjectName(QStringLiteral("main_vcfChrLbl"));
        main_vcfChrLbl->setEnabled(false);

        gridLayout->addWidget(main_vcfChrLbl, 4, 1, 1, 1);


        verticalLayout_7->addLayout(gridLayout);


        InputLayout1->addWidget(main_vcfGrp);

        main_sampleGrp = new QGroupBox(mainTab);
        main_sampleGrp->setObjectName(QStringLiteral("main_sampleGrp"));
        verticalLayout_2 = new QVBoxLayout(main_sampleGrp);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        SampleDirLayout = new QHBoxLayout();
        SampleDirLayout->setSpacing(6);
        SampleDirLayout->setObjectName(QStringLiteral("SampleDirLayout"));
        main_sampleDirTxt = new QLineEdit(main_sampleGrp);
        main_sampleDirTxt->setObjectName(QStringLiteral("main_sampleDirTxt"));
        main_sampleDirTxt->setReadOnly(true);

        SampleDirLayout->addWidget(main_sampleDirTxt);

        main_sampleDirBtn = new QPushButton(main_sampleGrp);
        main_sampleDirBtn->setObjectName(QStringLiteral("main_sampleDirBtn"));
        main_sampleDirBtn->setMaximumSize(QSize(20, 16777215));

        SampleDirLayout->addWidget(main_sampleDirBtn);


        verticalLayout_2->addLayout(SampleDirLayout);

        SampleParamLayout = new QHBoxLayout();
        SampleParamLayout->setSpacing(6);
        SampleParamLayout->setObjectName(QStringLiteral("SampleParamLayout"));
        main_sampleDepthTxt = new QLineEdit(main_sampleGrp);
        main_sampleDepthTxt->setObjectName(QStringLiteral("main_sampleDepthTxt"));

        SampleParamLayout->addWidget(main_sampleDepthTxt);

        main_sampleDepthLbl = new QLabel(main_sampleGrp);
        main_sampleDepthLbl->setObjectName(QStringLiteral("main_sampleDepthLbl"));

        SampleParamLayout->addWidget(main_sampleDepthLbl);


        verticalLayout_2->addLayout(SampleParamLayout);


        InputLayout1->addWidget(main_sampleGrp);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        InputLayout1->addItem(verticalSpacer_2);

        InputLayout1->setStretch(1, 15);

        horizontalLayout_2->addLayout(InputLayout1);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer);

        verticalLayout_5 = new QVBoxLayout();
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setObjectName(QStringLiteral("verticalLayout_5"));
        main_bedGrp = new QGroupBox(mainTab);
        main_bedGrp->setObjectName(QStringLiteral("main_bedGrp"));
        verticalLayout = new QVBoxLayout(main_bedGrp);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        BEDDirLayout = new QHBoxLayout();
        BEDDirLayout->setSpacing(6);
        BEDDirLayout->setObjectName(QStringLiteral("BEDDirLayout"));
        main_bedDirTxt = new QLineEdit(main_bedGrp);
        main_bedDirTxt->setObjectName(QStringLiteral("main_bedDirTxt"));
        main_bedDirTxt->setReadOnly(true);

        BEDDirLayout->addWidget(main_bedDirTxt);

        main_bedDirBtn = new QPushButton(main_bedGrp);
        main_bedDirBtn->setObjectName(QStringLiteral("main_bedDirBtn"));
        main_bedDirBtn->setMaximumSize(QSize(20, 16777215));

        BEDDirLayout->addWidget(main_bedDirBtn);


        verticalLayout->addLayout(BEDDirLayout);

        main_bedCollapseGrp = new QGroupBox(main_bedGrp);
        main_bedCollapseGrp->setObjectName(QStringLiteral("main_bedCollapseGrp"));
        horizontalLayout_13 = new QHBoxLayout(main_bedCollapseGrp);
        horizontalLayout_13->setSpacing(6);
        horizontalLayout_13->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_13->setObjectName(QStringLiteral("horizontalLayout_13"));
        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        main_bedCollapseExonBtn = new QRadioButton(main_bedCollapseGrp);
        main_bedCollapseExonBtn->setObjectName(QStringLiteral("main_bedCollapseExonBtn"));

        horizontalLayout_5->addWidget(main_bedCollapseExonBtn);

        main_bedCollapseGeneBtn = new QRadioButton(main_bedCollapseGrp);
        main_bedCollapseGeneBtn->setObjectName(QStringLiteral("main_bedCollapseGeneBtn"));
        main_bedCollapseGeneBtn->setChecked(false);

        horizontalLayout_5->addWidget(main_bedCollapseGeneBtn);

        main_bedCollapseKBtn = new QRadioButton(main_bedCollapseGrp);
        main_bedCollapseKBtn->setObjectName(QStringLiteral("main_bedCollapseKBtn"));
        main_bedCollapseKBtn->setChecked(true);

        horizontalLayout_5->addWidget(main_bedCollapseKBtn);

        main_bedCollapseKTxt = new QLineEdit(main_bedCollapseGrp);
        main_bedCollapseKTxt->setObjectName(QStringLiteral("main_bedCollapseKTxt"));

        horizontalLayout_5->addWidget(main_bedCollapseKTxt);


        horizontalLayout_13->addLayout(horizontalLayout_5);


        verticalLayout->addWidget(main_bedCollapseGrp);


        verticalLayout_5->addWidget(main_bedGrp);

        main_testGrp = new QGroupBox(mainTab);
        main_testGrp->setObjectName(QStringLiteral("main_testGrp"));
        horizontalLayout_12 = new QHBoxLayout(main_testGrp);
        horizontalLayout_12->setSpacing(6);
        horizontalLayout_12->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_12->setObjectName(QStringLiteral("horizontalLayout_12"));
        TestLayout = new QVBoxLayout();
        TestLayout->setSpacing(6);
        TestLayout->setObjectName(QStringLiteral("TestLayout"));
        main_testCommonBtn = new QRadioButton(main_testGrp);
        main_testCommonBtn->setObjectName(QStringLiteral("main_testCommonBtn"));
        main_testCommonBtn->setCheckable(true);
        main_testCommonBtn->setChecked(true);

        TestLayout->addWidget(main_testCommonBtn);

        main_testRareCastBtn = new QRadioButton(main_testGrp);
        main_testRareCastBtn->setObjectName(QStringLiteral("main_testRareCastBtn"));
        main_testRareCastBtn->setChecked(false);

        TestLayout->addWidget(main_testRareCastBtn);

        main_testRareCalphaBtn = new QRadioButton(main_testGrp);
        main_testRareCalphaBtn->setObjectName(QStringLiteral("main_testRareCalphaBtn"));

        TestLayout->addWidget(main_testRareCalphaBtn);


        horizontalLayout_12->addLayout(TestLayout);

        verticalLayout_4 = new QVBoxLayout();
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        main_testBootChk = new QCheckBox(main_testGrp);
        main_testBootChk->setObjectName(QStringLiteral("main_testBootChk"));
        main_testBootChk->setEnabled(false);
        main_testBootChk->setChecked(false);

        verticalLayout_4->addWidget(main_testBootChk);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        main_testBootTxt = new QLineEdit(main_testGrp);
        main_testBootTxt->setObjectName(QStringLiteral("main_testBootTxt"));
        main_testBootTxt->setEnabled(false);

        horizontalLayout_7->addWidget(main_testBootTxt);

        main_testBootLbl = new QLabel(main_testGrp);
        main_testBootLbl->setObjectName(QStringLiteral("main_testBootLbl"));

        horizontalLayout_7->addWidget(main_testBootLbl);


        verticalLayout_4->addLayout(horizontalLayout_7);

        main_testStopChk = new QCheckBox(main_testGrp);
        main_testStopChk->setObjectName(QStringLiteral("main_testStopChk"));
        main_testStopChk->setEnabled(false);
        main_testStopChk->setChecked(true);

        verticalLayout_4->addWidget(main_testStopChk);


        horizontalLayout_12->addLayout(verticalLayout_4);


        verticalLayout_5->addWidget(main_testGrp);

        main_startGrp = new QGroupBox(mainTab);
        main_startGrp->setObjectName(QStringLiteral("main_startGrp"));
        horizontalLayout = new QHBoxLayout(main_startGrp);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        formLayout = new QFormLayout();
        formLayout->setSpacing(6);
        formLayout->setObjectName(QStringLiteral("formLayout"));
        main_threadsTxt = new QLineEdit(main_startGrp);
        main_threadsTxt->setObjectName(QStringLiteral("main_threadsTxt"));
        main_threadsTxt->setEnabled(true);

        formLayout->setWidget(0, QFormLayout::LabelRole, main_threadsTxt);

        main_threadsLbl = new QLabel(main_startGrp);
        main_threadsLbl->setObjectName(QStringLiteral("main_threadsLbl"));

        formLayout->setWidget(0, QFormLayout::FieldRole, main_threadsLbl);

        main_batchSizeTxt = new QLineEdit(main_startGrp);
        main_batchSizeTxt->setObjectName(QStringLiteral("main_batchSizeTxt"));
        main_batchSizeTxt->setEnabled(true);

        formLayout->setWidget(1, QFormLayout::LabelRole, main_batchSizeTxt);

        main_batchSizeLbl = new QLabel(main_startGrp);
        main_batchSizeLbl->setObjectName(QStringLiteral("main_batchSizeLbl"));

        formLayout->setWidget(1, QFormLayout::FieldRole, main_batchSizeLbl);

        main_plotChk = new QCheckBox(main_startGrp);
        main_plotChk->setObjectName(QStringLiteral("main_plotChk"));
        main_plotChk->setChecked(true);

        formLayout->setWidget(2, QFormLayout::LabelRole, main_plotChk);

        main_rvsChk = new QCheckBox(main_startGrp);
        main_rvsChk->setObjectName(QStringLiteral("main_rvsChk"));
        main_rvsChk->setChecked(true);

        formLayout->setWidget(3, QFormLayout::LabelRole, main_rvsChk);

        main_gtChk = new QCheckBox(main_startGrp);
        main_gtChk->setObjectName(QStringLiteral("main_gtChk"));
        main_gtChk->setChecked(false);
        main_gtChk->setTristate(false);

        formLayout->setWidget(4, QFormLayout::LabelRole, main_gtChk);


        horizontalLayout->addLayout(formLayout);

        main_runBtn = new QPushButton(main_startGrp);
        main_runBtn->setObjectName(QStringLiteral("main_runBtn"));
        QSizePolicy sizePolicy1(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(1);
        sizePolicy1.setVerticalStretch(1);
        sizePolicy1.setHeightForWidth(main_runBtn->sizePolicy().hasHeightForWidth());
        main_runBtn->setSizePolicy(sizePolicy1);

        horizontalLayout->addWidget(main_runBtn);

        main_stopBtn = new QPushButton(main_startGrp);
        main_stopBtn->setObjectName(QStringLiteral("main_stopBtn"));

        horizontalLayout->addWidget(main_stopBtn);

        main_randomBtn = new QPushButton(main_startGrp);
        main_randomBtn->setObjectName(QStringLiteral("main_randomBtn"));

        horizontalLayout->addWidget(main_randomBtn);

        main_nrandomTxt = new QLineEdit(main_startGrp);
        main_nrandomTxt->setObjectName(QStringLiteral("main_nrandomTxt"));

        horizontalLayout->addWidget(main_nrandomTxt);


        verticalLayout_5->addWidget(main_startGrp);


        horizontalLayout_2->addLayout(verticalLayout_5);

        horizontalLayout_2->setStretch(0, 100);
        horizontalLayout_2->setStretch(1, 5);
        tabWidget->addTab(mainTab, QString());
        simulationTab = new QWidget();
        simulationTab->setObjectName(QStringLiteral("simulationTab"));
        horizontalLayout_8 = new QHBoxLayout(simulationTab);
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        sim_groupGrp = new QGroupBox(simulationTab);
        sim_groupGrp->setObjectName(QStringLiteral("sim_groupGrp"));
        verticalLayout_10 = new QVBoxLayout(sim_groupGrp);
        verticalLayout_10->setSpacing(6);
        verticalLayout_10->setContentsMargins(11, 11, 11, 11);
        verticalLayout_10->setObjectName(QStringLiteral("verticalLayout_10"));
        verticalLayout_11 = new QVBoxLayout();
        verticalLayout_11->setSpacing(6);
        verticalLayout_11->setObjectName(QStringLiteral("verticalLayout_11"));
        sim_groupAddBtn = new QPushButton(sim_groupGrp);
        sim_groupAddBtn->setObjectName(QStringLiteral("sim_groupAddBtn"));

        verticalLayout_11->addWidget(sim_groupAddBtn);

        sim_groupRemoveBtn = new QPushButton(sim_groupGrp);
        sim_groupRemoveBtn->setObjectName(QStringLiteral("sim_groupRemoveBtn"));

        verticalLayout_11->addWidget(sim_groupRemoveBtn);


        verticalLayout_10->addLayout(verticalLayout_11);

        sim_groupTbl = new QTableWidget(sim_groupGrp);
        if (sim_groupTbl->columnCount() < 5)
            sim_groupTbl->setColumnCount(5);
        QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
        sim_groupTbl->setHorizontalHeaderItem(0, __qtablewidgetitem);
        QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
        sim_groupTbl->setHorizontalHeaderItem(1, __qtablewidgetitem1);
        QTableWidgetItem *__qtablewidgetitem2 = new QTableWidgetItem();
        sim_groupTbl->setHorizontalHeaderItem(2, __qtablewidgetitem2);
        QTableWidgetItem *__qtablewidgetitem3 = new QTableWidgetItem();
        sim_groupTbl->setHorizontalHeaderItem(3, __qtablewidgetitem3);
        QTableWidgetItem *__qtablewidgetitem4 = new QTableWidgetItem();
        sim_groupTbl->setHorizontalHeaderItem(4, __qtablewidgetitem4);
        sim_groupTbl->setObjectName(QStringLiteral("sim_groupTbl"));
        sim_groupTbl->setEnabled(true);
        sim_groupTbl->horizontalHeader()->setStretchLastSection(true);
        sim_groupTbl->verticalHeader()->setStretchLastSection(false);

        verticalLayout_10->addWidget(sim_groupTbl);


        horizontalLayout_8->addWidget(sim_groupGrp);

        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        sim_variantGrp = new QGroupBox(simulationTab);
        sim_variantGrp->setObjectName(QStringLiteral("sim_variantGrp"));
        gridLayout_2 = new QGridLayout(sim_variantGrp);
        gridLayout_2->setSpacing(6);
        gridLayout_2->setContentsMargins(11, 11, 11, 11);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        sim_variantMafMaxTxt = new QLineEdit(sim_variantGrp);
        sim_variantMafMaxTxt->setObjectName(QStringLiteral("sim_variantMafMaxTxt"));

        gridLayout_2->addWidget(sim_variantMafMaxTxt, 2, 0, 1, 1);

        sim_variantSizeTxt = new QLineEdit(sim_variantGrp);
        sim_variantSizeTxt->setObjectName(QStringLiteral("sim_variantSizeTxt"));

        gridLayout_2->addWidget(sim_variantSizeTxt, 0, 0, 1, 1);

        sim_variantMafMinLbl = new QLabel(sim_variantGrp);
        sim_variantMafMinLbl->setObjectName(QStringLiteral("sim_variantMafMinLbl"));

        gridLayout_2->addWidget(sim_variantMafMinLbl, 1, 1, 1, 1);

        sim_variantSizeLbl = new QLabel(sim_variantGrp);
        sim_variantSizeLbl->setObjectName(QStringLiteral("sim_variantSizeLbl"));

        gridLayout_2->addWidget(sim_variantSizeLbl, 0, 1, 1, 1);

        sim_variantMafMaxLbl = new QLabel(sim_variantGrp);
        sim_variantMafMaxLbl->setObjectName(QStringLiteral("sim_variantMafMaxLbl"));

        gridLayout_2->addWidget(sim_variantMafMaxLbl, 2, 1, 1, 1);

        sim_variantMafMinTxt = new QLineEdit(sim_variantGrp);
        sim_variantMafMinTxt->setObjectName(QStringLiteral("sim_variantMafMinTxt"));

        gridLayout_2->addWidget(sim_variantMafMinTxt, 1, 0, 1, 1);


        verticalLayout_3->addWidget(sim_variantGrp);

        groupBox = new QGroupBox(simulationTab);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        gridLayout_6 = new QGridLayout(groupBox);
        gridLayout_6->setSpacing(6);
        gridLayout_6->setContentsMargins(11, 11, 11, 11);
        gridLayout_6->setObjectName(QStringLiteral("gridLayout_6"));
        sim_powerStepLbl = new QLabel(groupBox);
        sim_powerStepLbl->setObjectName(QStringLiteral("sim_powerStepLbl"));

        gridLayout_6->addWidget(sim_powerStepLbl, 1, 1, 1, 1);

        sim_groupHighLowLbl = new QLabel(groupBox);
        sim_groupHighLowLbl->setObjectName(QStringLiteral("sim_groupHighLowLbl"));

        gridLayout_6->addWidget(sim_groupHighLowLbl, 0, 1, 1, 1);

        sim_powerStepTxt = new QLineEdit(groupBox);
        sim_powerStepTxt->setObjectName(QStringLiteral("sim_powerStepTxt"));

        gridLayout_6->addWidget(sim_powerStepTxt, 1, 0, 1, 1);

        sim_groupHighLowTxt = new QLineEdit(groupBox);
        sim_groupHighLowTxt->setObjectName(QStringLiteral("sim_groupHighLowTxt"));

        gridLayout_6->addWidget(sim_groupHighLowTxt, 0, 0, 1, 1);

        sim_oddsRatioLbl = new QLabel(groupBox);
        sim_oddsRatioLbl->setObjectName(QStringLiteral("sim_oddsRatioLbl"));

        gridLayout_6->addWidget(sim_oddsRatioLbl, 2, 1, 1, 1);

        sim_oddsRatioTxt = new QLineEdit(groupBox);
        sim_oddsRatioTxt->setObjectName(QStringLiteral("sim_oddsRatioTxt"));

        gridLayout_6->addWidget(sim_oddsRatioTxt, 2, 0, 1, 1);


        verticalLayout_3->addWidget(groupBox);

        sim_testGrp = new QGroupBox(simulationTab);
        sim_testGrp->setObjectName(QStringLiteral("sim_testGrp"));
        horizontalLayout_4 = new QHBoxLayout(sim_testGrp);
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        verticalLayout_6 = new QVBoxLayout();
        verticalLayout_6->setSpacing(6);
        verticalLayout_6->setObjectName(QStringLiteral("verticalLayout_6"));
        sim_testCommonBtn = new QRadioButton(sim_testGrp);
        sim_testCommonBtn->setObjectName(QStringLiteral("sim_testCommonBtn"));
        sim_testCommonBtn->setChecked(true);

        verticalLayout_6->addWidget(sim_testCommonBtn);

        sim_testRareCastBtn = new QRadioButton(sim_testGrp);
        sim_testRareCastBtn->setObjectName(QStringLiteral("sim_testRareCastBtn"));
        sim_testRareCastBtn->setChecked(false);

        verticalLayout_6->addWidget(sim_testRareCastBtn);

        sim_testRareCalphaBtn = new QRadioButton(sim_testGrp);
        sim_testRareCalphaBtn->setObjectName(QStringLiteral("sim_testRareCalphaBtn"));

        verticalLayout_6->addWidget(sim_testRareCalphaBtn);


        horizontalLayout_4->addLayout(verticalLayout_6);

        verticalLayout_9 = new QVBoxLayout();
        verticalLayout_9->setSpacing(6);
        verticalLayout_9->setObjectName(QStringLiteral("verticalLayout_9"));
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        gridLayout_4 = new QGridLayout();
        gridLayout_4->setSpacing(6);
        gridLayout_4->setObjectName(QStringLiteral("gridLayout_4"));
        sim_testBootLbl = new QLabel(sim_testGrp);
        sim_testBootLbl->setObjectName(QStringLiteral("sim_testBootLbl"));
        sim_testBootLbl->setEnabled(false);

        gridLayout_4->addWidget(sim_testBootLbl, 0, 1, 1, 1);

        sim_collapseLbl = new QLabel(sim_testGrp);
        sim_collapseLbl->setObjectName(QStringLiteral("sim_collapseLbl"));
        sim_collapseLbl->setEnabled(false);

        gridLayout_4->addWidget(sim_collapseLbl, 1, 1, 1, 1);

        sim_testBootTxt = new QLineEdit(sim_testGrp);
        sim_testBootTxt->setObjectName(QStringLiteral("sim_testBootTxt"));
        sim_testBootTxt->setEnabled(false);

        gridLayout_4->addWidget(sim_testBootTxt, 0, 0, 1, 1);

        sim_collapseTxt = new QLineEdit(sim_testGrp);
        sim_collapseTxt->setObjectName(QStringLiteral("sim_collapseTxt"));
        sim_collapseTxt->setEnabled(false);

        gridLayout_4->addWidget(sim_collapseTxt, 1, 0, 1, 1);


        horizontalLayout_3->addLayout(gridLayout_4);


        verticalLayout_9->addLayout(horizontalLayout_3);

        sim_testStopChk = new QCheckBox(sim_testGrp);
        sim_testStopChk->setObjectName(QStringLiteral("sim_testStopChk"));
        sim_testStopChk->setEnabled(false);

        verticalLayout_9->addWidget(sim_testStopChk);


        horizontalLayout_4->addLayout(verticalLayout_9);


        verticalLayout_3->addWidget(sim_testGrp);

        sim_runBtn = new QPushButton(simulationTab);
        sim_runBtn->setObjectName(QStringLiteral("sim_runBtn"));
        sim_runBtn->setAutoDefault(false);
        sim_runBtn->setFlat(false);

        verticalLayout_3->addWidget(sim_runBtn);

        sim_stopBtn = new QPushButton(simulationTab);
        sim_stopBtn->setObjectName(QStringLiteral("sim_stopBtn"));
        sim_stopBtn->setEnabled(false);

        verticalLayout_3->addWidget(sim_stopBtn);


        horizontalLayout_8->addLayout(verticalLayout_3);

        horizontalLayout_8->setStretch(0, 4);
        horizontalLayout_8->setStretch(1, 1);
        tabWidget->addTab(simulationTab, QString());
        qSimTab = new QWidget();
        qSimTab->setObjectName(QStringLiteral("qSimTab"));
        horizontalLayout_10 = new QHBoxLayout(qSimTab);
        horizontalLayout_10->setSpacing(6);
        horizontalLayout_10->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_10->setObjectName(QStringLiteral("horizontalLayout_10"));
        qsim_groupGrp = new QGroupBox(qSimTab);
        qsim_groupGrp->setObjectName(QStringLiteral("qsim_groupGrp"));
        verticalLayout_12 = new QVBoxLayout(qsim_groupGrp);
        verticalLayout_12->setSpacing(6);
        verticalLayout_12->setContentsMargins(11, 11, 11, 11);
        verticalLayout_12->setObjectName(QStringLiteral("verticalLayout_12"));
        verticalLayout_13 = new QVBoxLayout();
        verticalLayout_13->setSpacing(6);
        verticalLayout_13->setObjectName(QStringLiteral("verticalLayout_13"));
        qsim_groupAddBtn = new QPushButton(qsim_groupGrp);
        qsim_groupAddBtn->setObjectName(QStringLiteral("qsim_groupAddBtn"));

        verticalLayout_13->addWidget(qsim_groupAddBtn);

        qsim_groupRemoveBtn = new QPushButton(qsim_groupGrp);
        qsim_groupRemoveBtn->setObjectName(QStringLiteral("qsim_groupRemoveBtn"));

        verticalLayout_13->addWidget(qsim_groupRemoveBtn);


        verticalLayout_12->addLayout(verticalLayout_13);

        qsim_groupTbl = new QTableWidget(qsim_groupGrp);
        if (qsim_groupTbl->columnCount() < 5)
            qsim_groupTbl->setColumnCount(5);
        QTableWidgetItem *__qtablewidgetitem5 = new QTableWidgetItem();
        qsim_groupTbl->setHorizontalHeaderItem(0, __qtablewidgetitem5);
        QTableWidgetItem *__qtablewidgetitem6 = new QTableWidgetItem();
        qsim_groupTbl->setHorizontalHeaderItem(1, __qtablewidgetitem6);
        QTableWidgetItem *__qtablewidgetitem7 = new QTableWidgetItem();
        qsim_groupTbl->setHorizontalHeaderItem(2, __qtablewidgetitem7);
        QTableWidgetItem *__qtablewidgetitem8 = new QTableWidgetItem();
        qsim_groupTbl->setHorizontalHeaderItem(3, __qtablewidgetitem8);
        QTableWidgetItem *__qtablewidgetitem9 = new QTableWidgetItem();
        qsim_groupTbl->setHorizontalHeaderItem(4, __qtablewidgetitem9);
        qsim_groupTbl->setObjectName(QStringLiteral("qsim_groupTbl"));
        qsim_groupTbl->setEnabled(true);
        qsim_groupTbl->horizontalHeader()->setStretchLastSection(true);
        qsim_groupTbl->verticalHeader()->setStretchLastSection(false);

        verticalLayout_12->addWidget(qsim_groupTbl);


        horizontalLayout_10->addWidget(qsim_groupGrp);

        verticalLayout_14 = new QVBoxLayout();
        verticalLayout_14->setSpacing(6);
        verticalLayout_14->setObjectName(QStringLiteral("verticalLayout_14"));
        qsim_variantGrp = new QGroupBox(qSimTab);
        qsim_variantGrp->setObjectName(QStringLiteral("qsim_variantGrp"));
        gridLayout_3 = new QGridLayout(qsim_variantGrp);
        gridLayout_3->setSpacing(6);
        gridLayout_3->setContentsMargins(11, 11, 11, 11);
        gridLayout_3->setObjectName(QStringLiteral("gridLayout_3"));
        qsim_variantMafMaxTxt = new QLineEdit(qsim_variantGrp);
        qsim_variantMafMaxTxt->setObjectName(QStringLiteral("qsim_variantMafMaxTxt"));

        gridLayout_3->addWidget(qsim_variantMafMaxTxt, 2, 0, 1, 1);

        qsim_variantSizeTxt = new QLineEdit(qsim_variantGrp);
        qsim_variantSizeTxt->setObjectName(QStringLiteral("qsim_variantSizeTxt"));

        gridLayout_3->addWidget(qsim_variantSizeTxt, 0, 0, 1, 1);

        qsim_variantMafMinLbl = new QLabel(qsim_variantGrp);
        qsim_variantMafMinLbl->setObjectName(QStringLiteral("qsim_variantMafMinLbl"));

        gridLayout_3->addWidget(qsim_variantMafMinLbl, 1, 1, 1, 1);

        qsim_variantSizeLbl = new QLabel(qsim_variantGrp);
        qsim_variantSizeLbl->setObjectName(QStringLiteral("qsim_variantSizeLbl"));

        gridLayout_3->addWidget(qsim_variantSizeLbl, 0, 1, 1, 1);

        qsim_variantMafMaxLbl = new QLabel(qsim_variantGrp);
        qsim_variantMafMaxLbl->setObjectName(QStringLiteral("qsim_variantMafMaxLbl"));

        gridLayout_3->addWidget(qsim_variantMafMaxLbl, 2, 1, 1, 1);

        qsim_variantMafMinTxt = new QLineEdit(qsim_variantGrp);
        qsim_variantMafMinTxt->setObjectName(QStringLiteral("qsim_variantMafMinTxt"));

        gridLayout_3->addWidget(qsim_variantMafMinTxt, 1, 0, 1, 1);


        verticalLayout_14->addWidget(qsim_variantGrp);

        qsim_runGrp = new QGroupBox(qSimTab);
        qsim_runGrp->setObjectName(QStringLiteral("qsim_runGrp"));
        gridLayout_7 = new QGridLayout(qsim_runGrp);
        gridLayout_7->setSpacing(6);
        gridLayout_7->setContentsMargins(11, 11, 11, 11);
        gridLayout_7->setObjectName(QStringLiteral("gridLayout_7"));
        qsim_powerStepLbl = new QLabel(qsim_runGrp);
        qsim_powerStepLbl->setObjectName(QStringLiteral("qsim_powerStepLbl"));

        gridLayout_7->addWidget(qsim_powerStepLbl, 1, 1, 1, 1);

        qsim_groupHighLowLbl = new QLabel(qsim_runGrp);
        qsim_groupHighLowLbl->setObjectName(QStringLiteral("qsim_groupHighLowLbl"));

        gridLayout_7->addWidget(qsim_groupHighLowLbl, 0, 1, 1, 1);

        qsim_powerStepTxt = new QLineEdit(qsim_runGrp);
        qsim_powerStepTxt->setObjectName(QStringLiteral("qsim_powerStepTxt"));

        gridLayout_7->addWidget(qsim_powerStepTxt, 1, 0, 1, 1);

        qsim_groupHighLowTxt = new QLineEdit(qsim_runGrp);
        qsim_groupHighLowTxt->setObjectName(QStringLiteral("qsim_groupHighLowTxt"));

        gridLayout_7->addWidget(qsim_groupHighLowTxt, 0, 0, 1, 1);

        qsim_R2Lbl = new QLabel(qsim_runGrp);
        qsim_R2Lbl->setObjectName(QStringLiteral("qsim_R2Lbl"));

        gridLayout_7->addWidget(qsim_R2Lbl, 2, 1, 1, 1);

        qsim_R2Txt = new QLineEdit(qsim_runGrp);
        qsim_R2Txt->setObjectName(QStringLiteral("qsim_R2Txt"));

        gridLayout_7->addWidget(qsim_R2Txt, 2, 0, 1, 1);


        verticalLayout_14->addWidget(qsim_runGrp);

        qsim_testGrp = new QGroupBox(qSimTab);
        qsim_testGrp->setObjectName(QStringLiteral("qsim_testGrp"));
        horizontalLayout_6 = new QHBoxLayout(qsim_testGrp);
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        verticalLayout_15 = new QVBoxLayout();
        verticalLayout_15->setSpacing(6);
        verticalLayout_15->setObjectName(QStringLiteral("verticalLayout_15"));
        qsim_testCommonBtn = new QRadioButton(qsim_testGrp);
        qsim_testCommonBtn->setObjectName(QStringLiteral("qsim_testCommonBtn"));
        qsim_testCommonBtn->setChecked(true);

        verticalLayout_15->addWidget(qsim_testCommonBtn);

        qsim_testRareCastBtn = new QRadioButton(qsim_testGrp);
        qsim_testRareCastBtn->setObjectName(QStringLiteral("qsim_testRareCastBtn"));
        qsim_testRareCastBtn->setChecked(false);

        verticalLayout_15->addWidget(qsim_testRareCastBtn);

        qsim_testRareCalphaBtn = new QRadioButton(qsim_testGrp);
        qsim_testRareCalphaBtn->setObjectName(QStringLiteral("qsim_testRareCalphaBtn"));

        verticalLayout_15->addWidget(qsim_testRareCalphaBtn);


        horizontalLayout_6->addLayout(verticalLayout_15);

        verticalLayout_16 = new QVBoxLayout();
        verticalLayout_16->setSpacing(6);
        verticalLayout_16->setObjectName(QStringLiteral("verticalLayout_16"));
        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        gridLayout_5 = new QGridLayout();
        gridLayout_5->setSpacing(6);
        gridLayout_5->setObjectName(QStringLiteral("gridLayout_5"));
        qsim_testBootLbl = new QLabel(qsim_testGrp);
        qsim_testBootLbl->setObjectName(QStringLiteral("qsim_testBootLbl"));
        qsim_testBootLbl->setEnabled(false);

        gridLayout_5->addWidget(qsim_testBootLbl, 0, 1, 1, 1);

        qsim_collapseLbl = new QLabel(qsim_testGrp);
        qsim_collapseLbl->setObjectName(QStringLiteral("qsim_collapseLbl"));
        qsim_collapseLbl->setEnabled(false);

        gridLayout_5->addWidget(qsim_collapseLbl, 1, 1, 1, 1);

        qsim_testBootTxt = new QLineEdit(qsim_testGrp);
        qsim_testBootTxt->setObjectName(QStringLiteral("qsim_testBootTxt"));
        qsim_testBootTxt->setEnabled(false);

        gridLayout_5->addWidget(qsim_testBootTxt, 0, 0, 1, 1);

        qsim_collapseTxt = new QLineEdit(qsim_testGrp);
        qsim_collapseTxt->setObjectName(QStringLiteral("qsim_collapseTxt"));
        qsim_collapseTxt->setEnabled(false);

        gridLayout_5->addWidget(qsim_collapseTxt, 1, 0, 1, 1);


        horizontalLayout_9->addLayout(gridLayout_5);


        verticalLayout_16->addLayout(horizontalLayout_9);

        qsim_testStopChk = new QCheckBox(qsim_testGrp);
        qsim_testStopChk->setObjectName(QStringLiteral("qsim_testStopChk"));
        qsim_testStopChk->setEnabled(false);

        verticalLayout_16->addWidget(qsim_testStopChk);


        horizontalLayout_6->addLayout(verticalLayout_16);


        verticalLayout_14->addWidget(qsim_testGrp);

        qsim_runBtn = new QPushButton(qSimTab);
        qsim_runBtn->setObjectName(QStringLiteral("qsim_runBtn"));
        qsim_runBtn->setAutoDefault(false);
        qsim_runBtn->setFlat(false);

        verticalLayout_14->addWidget(qsim_runBtn);

        qsim_stopBtn = new QPushButton(qSimTab);
        qsim_stopBtn->setObjectName(QStringLiteral("qsim_stopBtn"));
        qsim_stopBtn->setEnabled(false);

        verticalLayout_14->addWidget(qsim_stopBtn);


        horizontalLayout_10->addLayout(verticalLayout_14);

        horizontalLayout_10->setStretch(0, 4);
        horizontalLayout_10->setStretch(1, 1);
        tabWidget->addTab(qSimTab, QString());

        verticalLayout_8->addWidget(tabWidget);

        outputBox = new QTextEdit(MainWindowFrame);
        outputBox->setObjectName(QStringLiteral("outputBox"));
        QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(40);
        sizePolicy2.setHeightForWidth(outputBox->sizePolicy().hasHeightForWidth());
        outputBox->setSizePolicy(sizePolicy2);
        outputBox->setReadOnly(true);

        verticalLayout_8->addWidget(outputBox);

        MainWindow->setCentralWidget(MainWindowFrame);

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(1);
        sim_runBtn->setDefault(false);
        qsim_runBtn->setDefault(false);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "vikNGS GUI interface", nullptr));
        main_vcfGrp->setTitle(QApplication::translate("MainWindow", "Multisample VCF File", nullptr));
        main_vcfDirTxt->setText(QApplication::translate("MainWindow", "/media/scott/Rotom/chr2.merged.vcf", nullptr));
        main_vcfDirBtn->setText(QApplication::translate("MainWindow", "...", nullptr));
        main_vcfMissingLbl->setText(QApplication::translate("MainWindow", "Missing Threshold", nullptr));
        main_vcfWholeFileChk->setText(QApplication::translate("MainWindow", "Whole file", nullptr));
        main_vcfChromFilterTxt->setText(QString());
        main_vcfToPosLbl->setText(QApplication::translate("MainWindow", "End at", nullptr));
        main_vcfMissingTxt->setText(QApplication::translate("MainWindow", "0.1", nullptr));
#ifndef QT_NO_TOOLTIP
        main_vcfPassChk->setToolTip(QApplication::translate("MainWindow", "QUAL column in VCF must contain PASS", nullptr));
#endif // QT_NO_TOOLTIP
        main_vcfPassChk->setText(QApplication::translate("MainWindow", "Must PASS", nullptr));
        main_vcfFromPosLbl->setText(QApplication::translate("MainWindow", "Start from", nullptr));
        main_vcfMafLbl->setText(QApplication::translate("MainWindow", "MAF Cutoff", nullptr));
        main_vcfMafTxt->setText(QApplication::translate("MainWindow", "0.05", nullptr));
        main_vcfMafTxt->setPlaceholderText(QApplication::translate("MainWindow", "0.05", nullptr));
        main_vcfChrLbl->setText(QApplication::translate("MainWindow", "Chromosome", nullptr));
        main_sampleGrp->setTitle(QApplication::translate("MainWindow", "Sample Data File", nullptr));
        main_sampleDirTxt->setText(QApplication::translate("MainWindow", "/home/scott/vikNGS/vcf/chr22.txt", nullptr));
        main_sampleDirBtn->setText(QApplication::translate("MainWindow", "...", nullptr));
        main_sampleDepthTxt->setText(QApplication::translate("MainWindow", "30", nullptr));
        main_sampleDepthLbl->setText(QApplication::translate("MainWindow", "Read Depth High/Low Cutoff", nullptr));
        main_bedGrp->setTitle(QApplication::translate("MainWindow", "BED file (optional)", nullptr));
        main_bedDirTxt->setText(QString());
        main_bedDirBtn->setText(QApplication::translate("MainWindow", "...", nullptr));
        main_bedCollapseGrp->setTitle(QApplication::translate("MainWindow", "Collapse By", nullptr));
        main_bedCollapseExonBtn->setText(QApplication::translate("MainWindow", "Exon only", nullptr));
        main_bedCollapseGeneBtn->setText(QApplication::translate("MainWindow", "Entire gene", nullptr));
        main_bedCollapseKBtn->setText(QApplication::translate("MainWindow", "Every k ", nullptr));
        main_bedCollapseKTxt->setText(QApplication::translate("MainWindow", "5", nullptr));
        main_testGrp->setTitle(QApplication::translate("MainWindow", "Test", nullptr));
        main_testCommonBtn->setText(QApplication::translate("MainWindow", "Common", nullptr));
        main_testRareCastBtn->setText(QApplication::translate("MainWindow", "Rare (CAST)", nullptr));
        main_testRareCalphaBtn->setText(QApplication::translate("MainWindow", "Rare (c-alpha)", nullptr));
        main_testBootChk->setText(QApplication::translate("MainWindow", "Bootstrap", nullptr));
        main_testBootTxt->setText(QApplication::translate("MainWindow", "10000000", nullptr));
        main_testBootLbl->setText(QApplication::translate("MainWindow", "# iterations", nullptr));
        main_testStopChk->setText(QApplication::translate("MainWindow", "Stop Early", nullptr));
        main_startGrp->setTitle(QApplication::translate("MainWindow", "Start", nullptr));
        main_threadsTxt->setText(QApplication::translate("MainWindow", "2", nullptr));
        main_threadsLbl->setText(QApplication::translate("MainWindow", "# threads", nullptr));
        main_batchSizeTxt->setText(QApplication::translate("MainWindow", "100", nullptr));
        main_batchSizeLbl->setText(QApplication::translate("MainWindow", "Batch Size", nullptr));
        main_plotChk->setText(QApplication::translate("MainWindow", "Plot results", nullptr));
        main_rvsChk->setText(QApplication::translate("MainWindow", "RVS", nullptr));
        main_gtChk->setText(QApplication::translate("MainWindow", "Use GT", nullptr));
        main_runBtn->setText(QApplication::translate("MainWindow", "RUN", nullptr));
        main_stopBtn->setText(QApplication::translate("MainWindow", "STOP", nullptr));
        main_randomBtn->setText(QApplication::translate("MainWindow", "Random", nullptr));
        main_nrandomTxt->setText(QApplication::translate("MainWindow", "1000", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(mainTab), QApplication::translate("MainWindow", "Association Test", nullptr));
        sim_groupGrp->setTitle(QApplication::translate("MainWindow", "Group", nullptr));
        sim_groupAddBtn->setText(QApplication::translate("MainWindow", "Add Group", nullptr));
        sim_groupRemoveBtn->setText(QApplication::translate("MainWindow", "Remove Selected", nullptr));
        QTableWidgetItem *___qtablewidgetitem = sim_groupTbl->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QApplication::translate("MainWindow", "# Individuals", nullptr));
        QTableWidgetItem *___qtablewidgetitem1 = sim_groupTbl->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QApplication::translate("MainWindow", "Cohort", nullptr));
        QTableWidgetItem *___qtablewidgetitem2 = sim_groupTbl->horizontalHeaderItem(2);
        ___qtablewidgetitem2->setText(QApplication::translate("MainWindow", "Mean Depth", nullptr));
        QTableWidgetItem *___qtablewidgetitem3 = sim_groupTbl->horizontalHeaderItem(3);
        ___qtablewidgetitem3->setText(QApplication::translate("MainWindow", "Depth SD", nullptr));
        QTableWidgetItem *___qtablewidgetitem4 = sim_groupTbl->horizontalHeaderItem(4);
        ___qtablewidgetitem4->setText(QApplication::translate("MainWindow", "Error Rate", nullptr));
        sim_variantGrp->setTitle(QApplication::translate("MainWindow", "Variant Parameters", nullptr));
        sim_variantMafMaxTxt->setText(QApplication::translate("MainWindow", "0.45", nullptr));
        sim_variantSizeTxt->setText(QApplication::translate("MainWindow", "1000", nullptr));
#ifndef QT_NO_TOOLTIP
        sim_variantMafMinLbl->setToolTip(QApplication::translate("MainWindow", "Minor allele frequency", nullptr));
#endif // QT_NO_TOOLTIP
        sim_variantMafMinLbl->setText(QApplication::translate("MainWindow", "MAF min", nullptr));
#ifndef QT_NO_TOOLTIP
        sim_variantSizeLbl->setToolTip(QApplication::translate("MainWindow", "Number of variants to test", nullptr));
#endif // QT_NO_TOOLTIP
        sim_variantSizeLbl->setText(QApplication::translate("MainWindow", "# of variants", nullptr));
#ifndef QT_NO_TOOLTIP
        sim_variantMafMaxLbl->setToolTip(QApplication::translate("MainWindow", "Minor allele frequency", nullptr));
#endif // QT_NO_TOOLTIP
        sim_variantMafMaxLbl->setText(QApplication::translate("MainWindow", "MAF max", nullptr));
        sim_variantMafMinTxt->setText(QApplication::translate("MainWindow", "0.1", nullptr));
        groupBox->setTitle(QApplication::translate("MainWindow", "Run Parameters", nullptr));
#ifndef QT_NO_TOOLTIP
        sim_powerStepLbl->setToolTip(QApplication::translate("MainWindow", "Mean read depth of sequencing group", nullptr));
#endif // QT_NO_TOOLTIP
        sim_powerStepLbl->setText(QApplication::translate("MainWindow", "Steps", nullptr));
#ifndef QT_NO_TOOLTIP
        sim_groupHighLowLbl->setToolTip(QApplication::translate("MainWindow", "Any mean read depth less than this will be considered a low read group", nullptr));
#endif // QT_NO_TOOLTIP
        sim_groupHighLowLbl->setText(QApplication::translate("MainWindow", "High/low cutoff", nullptr));
        sim_powerStepTxt->setText(QApplication::translate("MainWindow", "5", nullptr));
        sim_groupHighLowTxt->setText(QApplication::translate("MainWindow", "30", nullptr));
#ifndef QT_NO_TOOLTIP
        sim_oddsRatioLbl->setToolTip(QApplication::translate("MainWindow", "Odds ratio of disease", nullptr));
#endif // QT_NO_TOOLTIP
        sim_oddsRatioLbl->setText(QApplication::translate("MainWindow", "Odds Ratio", nullptr));
        sim_oddsRatioTxt->setText(QApplication::translate("MainWindow", "1.0", nullptr));
        sim_testGrp->setTitle(QApplication::translate("MainWindow", "Test Parameters", nullptr));
        sim_testCommonBtn->setText(QApplication::translate("MainWindow", "Common", nullptr));
        sim_testRareCastBtn->setText(QApplication::translate("MainWindow", "Rare (CAST)", nullptr));
        sim_testRareCalphaBtn->setText(QApplication::translate("MainWindow", "Rare (c-alpha)", nullptr));
        sim_testBootLbl->setText(QApplication::translate("MainWindow", "# iterations", nullptr));
        sim_collapseLbl->setText(QApplication::translate("MainWindow", "Collapse ", nullptr));
        sim_testBootTxt->setText(QApplication::translate("MainWindow", "100", nullptr));
        sim_collapseTxt->setText(QApplication::translate("MainWindow", "5", nullptr));
        sim_testStopChk->setText(QApplication::translate("MainWindow", "Stop Early", nullptr));
        sim_runBtn->setText(QApplication::translate("MainWindow", "RUN", nullptr));
        sim_stopBtn->setText(QApplication::translate("MainWindow", "STOP", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(simulationTab), QApplication::translate("MainWindow", "Power Simulation (case/control)", nullptr));
        qsim_groupGrp->setTitle(QApplication::translate("MainWindow", "Group", nullptr));
        qsim_groupAddBtn->setText(QApplication::translate("MainWindow", "Add Group", nullptr));
        qsim_groupRemoveBtn->setText(QApplication::translate("MainWindow", "Remove Selected", nullptr));
        QTableWidgetItem *___qtablewidgetitem5 = qsim_groupTbl->horizontalHeaderItem(0);
        ___qtablewidgetitem5->setText(QApplication::translate("MainWindow", "# Individuals", nullptr));
        QTableWidgetItem *___qtablewidgetitem6 = qsim_groupTbl->horizontalHeaderItem(1);
        ___qtablewidgetitem6->setText(QApplication::translate("MainWindow", "Cohort", nullptr));
        QTableWidgetItem *___qtablewidgetitem7 = qsim_groupTbl->horizontalHeaderItem(2);
        ___qtablewidgetitem7->setText(QApplication::translate("MainWindow", "Mean Depth", nullptr));
        QTableWidgetItem *___qtablewidgetitem8 = qsim_groupTbl->horizontalHeaderItem(3);
        ___qtablewidgetitem8->setText(QApplication::translate("MainWindow", "Depth SD", nullptr));
        QTableWidgetItem *___qtablewidgetitem9 = qsim_groupTbl->horizontalHeaderItem(4);
        ___qtablewidgetitem9->setText(QApplication::translate("MainWindow", "Error Rate", nullptr));
        qsim_variantGrp->setTitle(QApplication::translate("MainWindow", "Variant Parameters", nullptr));
        qsim_variantMafMaxTxt->setText(QApplication::translate("MainWindow", "0.45", nullptr));
        qsim_variantSizeTxt->setText(QApplication::translate("MainWindow", "1000", nullptr));
#ifndef QT_NO_TOOLTIP
        qsim_variantMafMinLbl->setToolTip(QApplication::translate("MainWindow", "Minor allele frequency", nullptr));
#endif // QT_NO_TOOLTIP
        qsim_variantMafMinLbl->setText(QApplication::translate("MainWindow", "MAF min", nullptr));
#ifndef QT_NO_TOOLTIP
        qsim_variantSizeLbl->setToolTip(QApplication::translate("MainWindow", "Number of variants to test", nullptr));
#endif // QT_NO_TOOLTIP
        qsim_variantSizeLbl->setText(QApplication::translate("MainWindow", "# of variants", nullptr));
#ifndef QT_NO_TOOLTIP
        qsim_variantMafMaxLbl->setToolTip(QApplication::translate("MainWindow", "Minor allele frequency", nullptr));
#endif // QT_NO_TOOLTIP
        qsim_variantMafMaxLbl->setText(QApplication::translate("MainWindow", "MAF max", nullptr));
        qsim_variantMafMinTxt->setText(QApplication::translate("MainWindow", "0.1", nullptr));
        qsim_runGrp->setTitle(QApplication::translate("MainWindow", "Run Parameters", nullptr));
#ifndef QT_NO_TOOLTIP
        qsim_powerStepLbl->setToolTip(QApplication::translate("MainWindow", "Mean read depth of sequencing group", nullptr));
#endif // QT_NO_TOOLTIP
        qsim_powerStepLbl->setText(QApplication::translate("MainWindow", "Steps", nullptr));
#ifndef QT_NO_TOOLTIP
        qsim_groupHighLowLbl->setToolTip(QApplication::translate("MainWindow", "Any mean read depth less than this will be considered a low read group", nullptr));
#endif // QT_NO_TOOLTIP
        qsim_groupHighLowLbl->setText(QApplication::translate("MainWindow", "High/low cutoff", nullptr));
        qsim_powerStepTxt->setText(QApplication::translate("MainWindow", "5", nullptr));
        qsim_groupHighLowTxt->setText(QApplication::translate("MainWindow", "30", nullptr));
#ifndef QT_NO_TOOLTIP
        qsim_R2Lbl->setToolTip(QApplication::translate("MainWindow", "Odds ratio of disease", nullptr));
#endif // QT_NO_TOOLTIP
        qsim_R2Lbl->setText(QApplication::translate("MainWindow", "R2", nullptr));
        qsim_R2Txt->setText(QApplication::translate("MainWindow", "1.0", nullptr));
        qsim_testGrp->setTitle(QApplication::translate("MainWindow", "Test Parameters", nullptr));
        qsim_testCommonBtn->setText(QApplication::translate("MainWindow", "Common", nullptr));
        qsim_testRareCastBtn->setText(QApplication::translate("MainWindow", "Rare (CAST)", nullptr));
        qsim_testRareCalphaBtn->setText(QApplication::translate("MainWindow", "Rare (c-alpha)", nullptr));
        qsim_testBootLbl->setText(QApplication::translate("MainWindow", "# iterations", nullptr));
        qsim_collapseLbl->setText(QApplication::translate("MainWindow", "Collapse ", nullptr));
        qsim_testBootTxt->setText(QApplication::translate("MainWindow", "100", nullptr));
        qsim_collapseTxt->setText(QApplication::translate("MainWindow", "5", nullptr));
        qsim_testStopChk->setText(QApplication::translate("MainWindow", "Stop Early", nullptr));
        qsim_runBtn->setText(QApplication::translate("MainWindow", "RUN", nullptr));
        qsim_stopBtn->setText(QApplication::translate("MainWindow", "STOP", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(qSimTab), QApplication::translate("MainWindow", "Power Simulation (quantitative)", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
