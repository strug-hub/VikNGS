/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.11.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../gui/src/windows/mainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#include <QtCore/QVector>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.11.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_MainWindow_t {
    QByteArrayData data[62];
    char stringdata0[1080];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MainWindow_t qt_meta_stringdata_MainWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "MainWindow"
QT_MOC_LITERAL(1, 11, 12), // "sendPlotData"
QT_MOC_LITERAL(2, 24, 0), // ""
QT_MOC_LITERAL(3, 25, 15), // "QVector<double>"
QT_MOC_LITERAL(4, 41, 6), // "values"
QT_MOC_LITERAL(5, 48, 11), // "jobFinished"
QT_MOC_LITERAL(6, 60, 16), // "QVector<Variant>"
QT_MOC_LITERAL(7, 77, 8), // "variants"
QT_MOC_LITERAL(8, 86, 11), // "printOutput"
QT_MOC_LITERAL(9, 98, 6), // "string"
QT_MOC_LITERAL(10, 105, 1), // "c"
QT_MOC_LITERAL(11, 107, 10), // "greyOutput"
QT_MOC_LITERAL(12, 118, 15), // "constructGroups"
QT_MOC_LITERAL(13, 134, 35), // "std::vector<SimulationRequest..."
QT_MOC_LITERAL(14, 170, 6), // "ntests"
QT_MOC_LITERAL(15, 177, 16), // "constructRequest"
QT_MOC_LITERAL(16, 194, 17), // "SimulationRequest"
QT_MOC_LITERAL(17, 212, 6), // "groups"
QT_MOC_LITERAL(18, 219, 21), // "on_sim_runBtn_clicked"
QT_MOC_LITERAL(19, 241, 22), // "on_sim_stopBtn_clicked"
QT_MOC_LITERAL(20, 264, 8), // "addGroup"
QT_MOC_LITERAL(21, 273, 1), // "n"
QT_MOC_LITERAL(22, 275, 7), // "control"
QT_MOC_LITERAL(23, 283, 5), // "depth"
QT_MOC_LITERAL(24, 289, 7), // "sdDepth"
QT_MOC_LITERAL(25, 297, 9), // "errorRate"
QT_MOC_LITERAL(26, 307, 26), // "on_sim_groupAddBtn_clicked"
QT_MOC_LITERAL(27, 334, 29), // "on_sim_groupRemoveBtn_clicked"
QT_MOC_LITERAL(28, 364, 13), // "simEnableRare"
QT_MOC_LITERAL(29, 378, 5), // "valid"
QT_MOC_LITERAL(30, 384, 30), // "on_sim_testRareCastBtn_toggled"
QT_MOC_LITERAL(31, 415, 7), // "checked"
QT_MOC_LITERAL(32, 423, 32), // "on_sim_testRareCalphaBtn_toggled"
QT_MOC_LITERAL(33, 456, 18), // "simulationFinished"
QT_MOC_LITERAL(34, 475, 34), // "std::vector<std::vector<Varia..."
QT_MOC_LITERAL(35, 510, 4), // "reqs"
QT_MOC_LITERAL(36, 515, 7), // "stopJob"
QT_MOC_LITERAL(37, 523, 9), // "enableRun"
QT_MOC_LITERAL(38, 533, 10), // "disableRun"
QT_MOC_LITERAL(39, 544, 25), // "on_main_vcfDirBtn_clicked"
QT_MOC_LITERAL(40, 570, 28), // "on_main_sampleDirBtn_clicked"
QT_MOC_LITERAL(41, 599, 25), // "on_main_bedDirBtn_clicked"
QT_MOC_LITERAL(42, 625, 22), // "on_main_runBtn_clicked"
QT_MOC_LITERAL(43, 648, 31), // "on_main_testRareCastBtn_toggled"
QT_MOC_LITERAL(44, 680, 33), // "on_main_testRareCalphaBtn_tog..."
QT_MOC_LITERAL(45, 714, 32), // "on_main_testBootChk_stateChanged"
QT_MOC_LITERAL(46, 747, 4), // "arg1"
QT_MOC_LITERAL(47, 752, 27), // "on_main_testBootChk_toggled"
QT_MOC_LITERAL(48, 780, 31), // "on_main_vcfWholeFileChk_toggled"
QT_MOC_LITERAL(49, 812, 25), // "on_main_randomBtn_pressed"
QT_MOC_LITERAL(50, 838, 16), // "qConstructGroups"
QT_MOC_LITERAL(51, 855, 3), // "run"
QT_MOC_LITERAL(52, 859, 17), // "qConstructRequest"
QT_MOC_LITERAL(53, 877, 9), // "qAddGroup"
QT_MOC_LITERAL(54, 887, 14), // "qSimEnableRare"
QT_MOC_LITERAL(55, 902, 5), // "value"
QT_MOC_LITERAL(56, 908, 22), // "on_qsim_runBtn_clicked"
QT_MOC_LITERAL(57, 931, 27), // "on_qsim_groupAddBtn_clicked"
QT_MOC_LITERAL(58, 959, 30), // "on_qsim_groupRemoveBtn_clicked"
QT_MOC_LITERAL(59, 990, 23), // "on_qsim_stopBtn_clicked"
QT_MOC_LITERAL(60, 1014, 31), // "on_qsim_testRareCastBtn_toggled"
QT_MOC_LITERAL(61, 1046, 33) // "on_qsim_testRareCalphaBtn_tog..."

    },
    "MainWindow\0sendPlotData\0\0QVector<double>\0"
    "values\0jobFinished\0QVector<Variant>\0"
    "variants\0printOutput\0string\0c\0greyOutput\0"
    "constructGroups\0std::vector<SimulationRequestGroup>\0"
    "ntests\0constructRequest\0SimulationRequest\0"
    "groups\0on_sim_runBtn_clicked\0"
    "on_sim_stopBtn_clicked\0addGroup\0n\0"
    "control\0depth\0sdDepth\0errorRate\0"
    "on_sim_groupAddBtn_clicked\0"
    "on_sim_groupRemoveBtn_clicked\0"
    "simEnableRare\0valid\0on_sim_testRareCastBtn_toggled\0"
    "checked\0on_sim_testRareCalphaBtn_toggled\0"
    "simulationFinished\0"
    "std::vector<std::vector<Variant> >\0"
    "reqs\0stopJob\0enableRun\0disableRun\0"
    "on_main_vcfDirBtn_clicked\0"
    "on_main_sampleDirBtn_clicked\0"
    "on_main_bedDirBtn_clicked\0"
    "on_main_runBtn_clicked\0"
    "on_main_testRareCastBtn_toggled\0"
    "on_main_testRareCalphaBtn_toggled\0"
    "on_main_testBootChk_stateChanged\0arg1\0"
    "on_main_testBootChk_toggled\0"
    "on_main_vcfWholeFileChk_toggled\0"
    "on_main_randomBtn_pressed\0qConstructGroups\0"
    "run\0qConstructRequest\0qAddGroup\0"
    "qSimEnableRare\0value\0on_qsim_runBtn_clicked\0"
    "on_qsim_groupAddBtn_clicked\0"
    "on_qsim_groupRemoveBtn_clicked\0"
    "on_qsim_stopBtn_clicked\0"
    "on_qsim_testRareCastBtn_toggled\0"
    "on_qsim_testRareCalphaBtn_toggled"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      38,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,  204,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       5,    1,  207,    2, 0x0a /* Public */,
       8,    2,  210,    2, 0x0a /* Public */,
      11,    0,  215,    2, 0x0a /* Public */,
      12,    1,  216,    2, 0x08 /* Private */,
      15,    1,  219,    2, 0x08 /* Private */,
      18,    0,  222,    2, 0x08 /* Private */,
      19,    0,  223,    2, 0x08 /* Private */,
      20,    5,  224,    2, 0x08 /* Private */,
      26,    0,  235,    2, 0x08 /* Private */,
      27,    0,  236,    2, 0x08 /* Private */,
      28,    1,  237,    2, 0x08 /* Private */,
      30,    1,  240,    2, 0x08 /* Private */,
      32,    1,  243,    2, 0x08 /* Private */,
      33,    2,  246,    2, 0x08 /* Private */,
      36,    0,  251,    2, 0x08 /* Private */,
      37,    0,  252,    2, 0x08 /* Private */,
      38,    0,  253,    2, 0x08 /* Private */,
      39,    0,  254,    2, 0x08 /* Private */,
      40,    0,  255,    2, 0x08 /* Private */,
      41,    0,  256,    2, 0x08 /* Private */,
      42,    0,  257,    2, 0x08 /* Private */,
      43,    1,  258,    2, 0x08 /* Private */,
      44,    1,  261,    2, 0x08 /* Private */,
      45,    1,  264,    2, 0x08 /* Private */,
      47,    1,  267,    2, 0x08 /* Private */,
      48,    1,  270,    2, 0x08 /* Private */,
      49,    0,  273,    2, 0x08 /* Private */,
      50,    2,  274,    2, 0x08 /* Private */,
      52,    1,  279,    2, 0x08 /* Private */,
      53,    5,  282,    2, 0x08 /* Private */,
      54,    1,  293,    2, 0x08 /* Private */,
      56,    0,  296,    2, 0x08 /* Private */,
      57,    0,  297,    2, 0x08 /* Private */,
      58,    0,  298,    2, 0x08 /* Private */,
      59,    0,  299,    2, 0x08 /* Private */,
      60,    1,  300,    2, 0x08 /* Private */,
      61,    1,  303,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    4,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 6,    7,
    QMetaType::Void, QMetaType::QString, QMetaType::QColor,    9,   10,
    QMetaType::Void,
    0x80000000 | 13, QMetaType::Int,   14,
    0x80000000 | 16, 0x80000000 | 13,   17,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString, QMetaType::Bool, QMetaType::QString, QMetaType::QString, QMetaType::QString,   21,   22,   23,   24,   25,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,   29,
    QMetaType::Void, QMetaType::Bool,   31,
    QMetaType::Void, QMetaType::Bool,   31,
    QMetaType::Void, 0x80000000 | 34, 0x80000000 | 16,    7,   35,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,   31,
    QMetaType::Void, QMetaType::Bool,   31,
    QMetaType::Void, QMetaType::Int,   46,
    QMetaType::Void, QMetaType::Bool,   31,
    QMetaType::Void, QMetaType::Bool,   31,
    QMetaType::Void,
    0x80000000 | 13, QMetaType::Int, QMetaType::Int,   51,   14,
    0x80000000 | 16, 0x80000000 | 13,   17,
    QMetaType::Void, QMetaType::QString, QMetaType::Bool, QMetaType::QString, QMetaType::QString, QMetaType::QString,   21,   22,   23,   24,   25,
    QMetaType::Void, QMetaType::Bool,   55,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,   31,
    QMetaType::Void, QMetaType::Bool,   31,

       0        // eod
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MainWindow *_t = static_cast<MainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->sendPlotData((*reinterpret_cast< QVector<double>(*)>(_a[1]))); break;
        case 1: _t->jobFinished((*reinterpret_cast< QVector<Variant>(*)>(_a[1]))); break;
        case 2: _t->printOutput((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< QColor(*)>(_a[2]))); break;
        case 3: _t->greyOutput(); break;
        case 4: { std::vector<SimulationRequestGroup> _r = _t->constructGroups((*reinterpret_cast< int(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< std::vector<SimulationRequestGroup>*>(_a[0]) = std::move(_r); }  break;
        case 5: { SimulationRequest _r = _t->constructRequest((*reinterpret_cast< std::vector<SimulationRequestGroup>(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< SimulationRequest*>(_a[0]) = std::move(_r); }  break;
        case 6: _t->on_sim_runBtn_clicked(); break;
        case 7: _t->on_sim_stopBtn_clicked(); break;
        case 8: _t->addGroup((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2])),(*reinterpret_cast< QString(*)>(_a[3])),(*reinterpret_cast< QString(*)>(_a[4])),(*reinterpret_cast< QString(*)>(_a[5]))); break;
        case 9: _t->on_sim_groupAddBtn_clicked(); break;
        case 10: _t->on_sim_groupRemoveBtn_clicked(); break;
        case 11: _t->simEnableRare((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 12: _t->on_sim_testRareCastBtn_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 13: _t->on_sim_testRareCalphaBtn_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 14: _t->simulationFinished((*reinterpret_cast< std::vector<std::vector<Variant> >(*)>(_a[1])),(*reinterpret_cast< SimulationRequest(*)>(_a[2]))); break;
        case 15: _t->stopJob(); break;
        case 16: _t->enableRun(); break;
        case 17: _t->disableRun(); break;
        case 18: _t->on_main_vcfDirBtn_clicked(); break;
        case 19: _t->on_main_sampleDirBtn_clicked(); break;
        case 20: _t->on_main_bedDirBtn_clicked(); break;
        case 21: _t->on_main_runBtn_clicked(); break;
        case 22: _t->on_main_testRareCastBtn_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 23: _t->on_main_testRareCalphaBtn_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 24: _t->on_main_testBootChk_stateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 25: _t->on_main_testBootChk_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 26: _t->on_main_vcfWholeFileChk_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 27: _t->on_main_randomBtn_pressed(); break;
        case 28: { std::vector<SimulationRequestGroup> _r = _t->qConstructGroups((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])));
            if (_a[0]) *reinterpret_cast< std::vector<SimulationRequestGroup>*>(_a[0]) = std::move(_r); }  break;
        case 29: { SimulationRequest _r = _t->qConstructRequest((*reinterpret_cast< std::vector<SimulationRequestGroup>(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< SimulationRequest*>(_a[0]) = std::move(_r); }  break;
        case 30: _t->qAddGroup((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2])),(*reinterpret_cast< QString(*)>(_a[3])),(*reinterpret_cast< QString(*)>(_a[4])),(*reinterpret_cast< QString(*)>(_a[5]))); break;
        case 31: _t->qSimEnableRare((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 32: _t->on_qsim_runBtn_clicked(); break;
        case 33: _t->on_qsim_groupAddBtn_clicked(); break;
        case 34: _t->on_qsim_groupRemoveBtn_clicked(); break;
        case 35: _t->on_qsim_stopBtn_clicked(); break;
        case 36: _t->on_qsim_testRareCastBtn_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 37: _t->on_qsim_testRareCalphaBtn_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 0:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QVector<double> >(); break;
            }
            break;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (MainWindow::*)(QVector<double> );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MainWindow::sendPlotData)) {
                *result = 0;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow.data,
      qt_meta_data_MainWindow,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 38)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 38;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 38)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 38;
    }
    return _id;
}

// SIGNAL 0
void MainWindow::sendPlotData(QVector<double> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
