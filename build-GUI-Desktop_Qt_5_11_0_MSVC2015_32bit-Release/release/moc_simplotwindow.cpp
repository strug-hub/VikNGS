/****************************************************************************
** Meta object code from reading C++ file 'simplotwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.11.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../gui/src/windows/simplotwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'simplotwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.11.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_SimPlotWindow_t {
    QByteArrayData data[27];
    char stringdata0[389];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_SimPlotWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_SimPlotWindow_t qt_meta_stringdata_SimPlotWindow = {
    {
QT_MOC_LITERAL(0, 0, 13), // "SimPlotWindow"
QT_MOC_LITERAL(1, 14, 10), // "initialize"
QT_MOC_LITERAL(2, 25, 0), // ""
QT_MOC_LITERAL(3, 26, 35), // "std::vector<std::vector<Varia..."
QT_MOC_LITERAL(4, 62, 8), // "variants"
QT_MOC_LITERAL(5, 71, 18), // "SimulationRequest&"
QT_MOC_LITERAL(6, 90, 3), // "req"
QT_MOC_LITERAL(7, 94, 5), // "title"
QT_MOC_LITERAL(8, 100, 10), // "getPvalues"
QT_MOC_LITERAL(9, 111, 15), // "filterCollapsed"
QT_MOC_LITERAL(10, 127, 34), // "std::vector<std::vector<Varia..."
QT_MOC_LITERAL(11, 162, 1), // "k"
QT_MOC_LITERAL(12, 164, 9), // "buildPlot"
QT_MOC_LITERAL(13, 174, 11), // "buildLegend"
QT_MOC_LITERAL(14, 186, 16), // "updateSampleSize"
QT_MOC_LITERAL(15, 203, 5), // "index"
QT_MOC_LITERAL(16, 209, 14), // "mouseMovePlot1"
QT_MOC_LITERAL(17, 224, 12), // "QMouseEvent*"
QT_MOC_LITERAL(18, 237, 5), // "event"
QT_MOC_LITERAL(19, 243, 14), // "mouseMovePlot2"
QT_MOC_LITERAL(20, 258, 15), // "mouseClickPlot1"
QT_MOC_LITERAL(21, 274, 15), // "mouseClickPlot2"
QT_MOC_LITERAL(22, 290, 33), // "on_simplot_alphaDial_valueCha..."
QT_MOC_LITERAL(23, 324, 5), // "value"
QT_MOC_LITERAL(24, 330, 31), // "on_simplot_alphaTxt_textChanged"
QT_MOC_LITERAL(25, 362, 4), // "arg1"
QT_MOC_LITERAL(26, 367, 21) // "on_pushButton_pressed"

    },
    "SimPlotWindow\0initialize\0\0"
    "std::vector<std::vector<Variant> >&\0"
    "variants\0SimulationRequest&\0req\0title\0"
    "getPvalues\0filterCollapsed\0"
    "std::vector<std::vector<Variant> >\0k\0"
    "buildPlot\0buildLegend\0updateSampleSize\0"
    "index\0mouseMovePlot1\0QMouseEvent*\0"
    "event\0mouseMovePlot2\0mouseClickPlot1\0"
    "mouseClickPlot2\0on_simplot_alphaDial_valueChanged\0"
    "value\0on_simplot_alphaTxt_textChanged\0"
    "arg1\0on_pushButton_pressed"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_SimPlotWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      13,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    3,   79,    2, 0x0a /* Public */,
       8,    1,   86,    2, 0x0a /* Public */,
       9,    2,   89,    2, 0x0a /* Public */,
      12,    0,   94,    2, 0x0a /* Public */,
      13,    0,   95,    2, 0x0a /* Public */,
      14,    1,   96,    2, 0x0a /* Public */,
      16,    1,   99,    2, 0x0a /* Public */,
      19,    1,  102,    2, 0x0a /* Public */,
      20,    1,  105,    2, 0x0a /* Public */,
      21,    1,  108,    2, 0x0a /* Public */,
      22,    1,  111,    2, 0x08 /* Private */,
      24,    1,  114,    2, 0x08 /* Private */,
      26,    0,  117,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3, 0x80000000 | 5, QMetaType::QString,    4,    6,    7,
    QMetaType::Void, 0x80000000 | 3,    4,
    0x80000000 | 10, 0x80000000 | 3, QMetaType::Int,    4,   11,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   15,
    QMetaType::Void, 0x80000000 | 17,   18,
    QMetaType::Void, 0x80000000 | 17,   18,
    QMetaType::Void, 0x80000000 | 17,   18,
    QMetaType::Void, 0x80000000 | 17,   18,
    QMetaType::Void, QMetaType::Int,   23,
    QMetaType::Void, QMetaType::QString,   25,
    QMetaType::Void,

       0        // eod
};

void SimPlotWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        SimPlotWindow *_t = static_cast<SimPlotWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->initialize((*reinterpret_cast< std::vector<std::vector<Variant> >(*)>(_a[1])),(*reinterpret_cast< SimulationRequest(*)>(_a[2])),(*reinterpret_cast< QString(*)>(_a[3]))); break;
        case 1: _t->getPvalues((*reinterpret_cast< std::vector<std::vector<Variant> >(*)>(_a[1]))); break;
        case 2: { std::vector<std::vector<Variant> > _r = _t->filterCollapsed((*reinterpret_cast< std::vector<std::vector<Variant> >(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])));
            if (_a[0]) *reinterpret_cast< std::vector<std::vector<Variant> >*>(_a[0]) = std::move(_r); }  break;
        case 3: _t->buildPlot(); break;
        case 4: _t->buildLegend(); break;
        case 5: _t->updateSampleSize((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: _t->mouseMovePlot1((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 7: _t->mouseMovePlot2((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 8: _t->mouseClickPlot1((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 9: _t->mouseClickPlot2((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 10: _t->on_simplot_alphaDial_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: _t->on_simplot_alphaTxt_textChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 12: _t->on_pushButton_pressed(); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject SimPlotWindow::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_SimPlotWindow.data,
      qt_meta_data_SimPlotWindow,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *SimPlotWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *SimPlotWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_SimPlotWindow.stringdata0))
        return static_cast<void*>(this);
    return QWidget::qt_metacast(_clname);
}

int SimPlotWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 13)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 13;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 13)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 13;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
