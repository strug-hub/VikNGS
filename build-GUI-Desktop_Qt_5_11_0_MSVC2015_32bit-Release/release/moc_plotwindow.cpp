/****************************************************************************
** Meta object code from reading C++ file 'plotwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.11.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../gui/src/windows/plotwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#include <QtCore/QVector>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'plotwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.11.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_PlotWindow_t {
    QByteArrayData data[48];
    char stringdata0[573];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_PlotWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_PlotWindow_t qt_meta_stringdata_PlotWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "PlotWindow"
QT_MOC_LITERAL(1, 11, 10), // "initialize"
QT_MOC_LITERAL(2, 22, 0), // ""
QT_MOC_LITERAL(3, 23, 16), // "QVector<Variant>"
QT_MOC_LITERAL(4, 40, 6), // "values"
QT_MOC_LITERAL(5, 47, 5), // "title"
QT_MOC_LITERAL(6, 53, 1), // "n"
QT_MOC_LITERAL(7, 55, 17), // "createChromosomes"
QT_MOC_LITERAL(8, 73, 8), // "variants"
QT_MOC_LITERAL(9, 82, 20), // "setRandomChromosomes"
QT_MOC_LITERAL(10, 103, 24), // "generateRandomChromosome"
QT_MOC_LITERAL(11, 128, 10), // "Chromosome"
QT_MOC_LITERAL(12, 139, 11), // "std::string"
QT_MOC_LITERAL(13, 151, 5), // "chrom"
QT_MOC_LITERAL(14, 157, 6), // "maxPos"
QT_MOC_LITERAL(15, 164, 15), // "mouseMoveWindow"
QT_MOC_LITERAL(16, 180, 12), // "QMouseEvent*"
QT_MOC_LITERAL(17, 193, 5), // "event"
QT_MOC_LITERAL(18, 199, 15), // "mouseMoveGenome"
QT_MOC_LITERAL(19, 215, 16), // "mouseClickGenome"
QT_MOC_LITERAL(20, 232, 21), // "updateChromosomeRange"
QT_MOC_LITERAL(21, 254, 8), // "QCPRange"
QT_MOC_LITERAL(22, 263, 8), // "newRange"
QT_MOC_LITERAL(23, 272, 19), // "getChromUnderCursor"
QT_MOC_LITERAL(24, 292, 10), // "resetColor"
QT_MOC_LITERAL(25, 303, 7), // "chrName"
QT_MOC_LITERAL(26, 311, 19), // "mouseMoveChromosome"
QT_MOC_LITERAL(27, 331, 20), // "mouseClickChromosome"
QT_MOC_LITERAL(28, 352, 15), // "buildGenomePlot"
QT_MOC_LITERAL(29, 368, 19), // "buildChromosomePlot"
QT_MOC_LITERAL(30, 388, 27), // "updateVariantHighlightLayer"
QT_MOC_LITERAL(31, 416, 7), // "Variant"
QT_MOC_LITERAL(32, 424, 13), // "moveRectangle"
QT_MOC_LITERAL(33, 438, 12), // "QCPItemRect*"
QT_MOC_LITERAL(34, 451, 4), // "rect"
QT_MOC_LITERAL(35, 456, 5), // "lower"
QT_MOC_LITERAL(36, 462, 5), // "upper"
QT_MOC_LITERAL(37, 468, 17), // "updateVariantInfo"
QT_MOC_LITERAL(38, 486, 7), // "variant"
QT_MOC_LITERAL(39, 494, 14), // "getGraphByName"
QT_MOC_LITERAL(40, 509, 9), // "QCPGraph*"
QT_MOC_LITERAL(41, 519, 12), // "QCustomPlot*"
QT_MOC_LITERAL(42, 532, 4), // "plot"
QT_MOC_LITERAL(43, 537, 4), // "name"
QT_MOC_LITERAL(44, 542, 18), // "findClosestVariant"
QT_MOC_LITERAL(45, 561, 1), // "x"
QT_MOC_LITERAL(46, 563, 1), // "y"
QT_MOC_LITERAL(47, 565, 7) // "maxDist"

    },
    "PlotWindow\0initialize\0\0QVector<Variant>\0"
    "values\0title\0n\0createChromosomes\0"
    "variants\0setRandomChromosomes\0"
    "generateRandomChromosome\0Chromosome\0"
    "std::string\0chrom\0maxPos\0mouseMoveWindow\0"
    "QMouseEvent*\0event\0mouseMoveGenome\0"
    "mouseClickGenome\0updateChromosomeRange\0"
    "QCPRange\0newRange\0getChromUnderCursor\0"
    "resetColor\0chrName\0mouseMoveChromosome\0"
    "mouseClickChromosome\0buildGenomePlot\0"
    "buildChromosomePlot\0updateVariantHighlightLayer\0"
    "Variant\0moveRectangle\0QCPItemRect*\0"
    "rect\0lower\0upper\0updateVariantInfo\0"
    "variant\0getGraphByName\0QCPGraph*\0"
    "QCustomPlot*\0plot\0name\0findClosestVariant\0"
    "x\0y\0maxDist"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_PlotWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      22,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    2,  124,    2, 0x0a /* Public */,
       1,    2,  129,    2, 0x0a /* Public */,
       7,    1,  134,    2, 0x0a /* Public */,
       9,    1,  137,    2, 0x0a /* Public */,
      10,    3,  140,    2, 0x0a /* Public */,
      15,    1,  147,    2, 0x0a /* Public */,
      18,    1,  150,    2, 0x0a /* Public */,
      19,    1,  153,    2, 0x0a /* Public */,
      20,    1,  156,    2, 0x0a /* Public */,
      23,    1,  159,    2, 0x0a /* Public */,
      24,    1,  162,    2, 0x0a /* Public */,
      26,    1,  165,    2, 0x0a /* Public */,
      27,    1,  168,    2, 0x0a /* Public */,
      28,    0,  171,    2, 0x0a /* Public */,
      29,    1,  172,    2, 0x0a /* Public */,
      30,    1,  175,    2, 0x0a /* Public */,
      32,    4,  178,    2, 0x0a /* Public */,
      32,    3,  187,    2, 0x2a /* Public | MethodCloned */,
      32,    2,  194,    2, 0x2a /* Public | MethodCloned */,
      37,    1,  199,    2, 0x0a /* Public */,
      39,    2,  202,    2, 0x0a /* Public */,
      44,    3,  207,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3, QMetaType::QString,    4,    5,
    QMetaType::Void, QMetaType::Int, QMetaType::QString,    6,    5,
    QMetaType::Void, 0x80000000 | 3,    8,
    QMetaType::Void, QMetaType::Int,    6,
    0x80000000 | 11, QMetaType::Int, 0x80000000 | 12, QMetaType::Int,    6,   13,   14,
    QMetaType::Void, 0x80000000 | 16,   17,
    QMetaType::Void, 0x80000000 | 16,   17,
    QMetaType::Void, 0x80000000 | 16,   17,
    QMetaType::Void, 0x80000000 | 21,   22,
    QMetaType::QString, 0x80000000 | 16,   17,
    QMetaType::Void, QMetaType::QString,   25,
    QMetaType::Void, 0x80000000 | 16,   17,
    QMetaType::Void, 0x80000000 | 16,   17,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,   25,
    QMetaType::Void, 0x80000000 | 31,    8,
    QMetaType::Void, 0x80000000 | 33, QMetaType::QString, QMetaType::Double, QMetaType::Double,   34,   25,   35,   36,
    QMetaType::Void, 0x80000000 | 33, QMetaType::QString, QMetaType::Double,   34,   25,   35,
    QMetaType::Void, 0x80000000 | 33, QMetaType::QString,   34,   25,
    QMetaType::Void, 0x80000000 | 31,   38,
    0x80000000 | 40, 0x80000000 | 41, QMetaType::QString,   42,   43,
    0x80000000 | 31, QMetaType::Double, QMetaType::Double, QMetaType::Double,   45,   46,   47,

       0        // eod
};

void PlotWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        PlotWindow *_t = static_cast<PlotWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->initialize((*reinterpret_cast< QVector<Variant>(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
        case 1: _t->initialize((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
        case 2: _t->createChromosomes((*reinterpret_cast< QVector<Variant>(*)>(_a[1]))); break;
        case 3: _t->setRandomChromosomes((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: { Chromosome _r = _t->generateRandomChromosome((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< std::string(*)>(_a[2])),(*reinterpret_cast< int(*)>(_a[3])));
            if (_a[0]) *reinterpret_cast< Chromosome*>(_a[0]) = std::move(_r); }  break;
        case 5: _t->mouseMoveWindow((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 6: _t->mouseMoveGenome((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 7: _t->mouseClickGenome((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 8: _t->updateChromosomeRange((*reinterpret_cast< const QCPRange(*)>(_a[1]))); break;
        case 9: { QString _r = _t->getChromUnderCursor((*reinterpret_cast< QMouseEvent*(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< QString*>(_a[0]) = std::move(_r); }  break;
        case 10: _t->resetColor((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 11: _t->mouseMoveChromosome((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 12: _t->mouseClickChromosome((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 13: _t->buildGenomePlot(); break;
        case 14: _t->buildChromosomePlot((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 15: _t->updateVariantHighlightLayer((*reinterpret_cast< Variant(*)>(_a[1]))); break;
        case 16: _t->moveRectangle((*reinterpret_cast< QCPItemRect*(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])),(*reinterpret_cast< double(*)>(_a[4]))); break;
        case 17: _t->moveRectangle((*reinterpret_cast< QCPItemRect*(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 18: _t->moveRectangle((*reinterpret_cast< QCPItemRect*(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
        case 19: _t->updateVariantInfo((*reinterpret_cast< Variant(*)>(_a[1]))); break;
        case 20: { QCPGraph* _r = _t->getGraphByName((*reinterpret_cast< QCustomPlot*(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2])));
            if (_a[0]) *reinterpret_cast< QCPGraph**>(_a[0]) = std::move(_r); }  break;
        case 21: { Variant _r = _t->findClosestVariant((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])));
            if (_a[0]) *reinterpret_cast< Variant*>(_a[0]) = std::move(_r); }  break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 16:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QCPItemRect* >(); break;
            }
            break;
        case 17:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QCPItemRect* >(); break;
            }
            break;
        case 18:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QCPItemRect* >(); break;
            }
            break;
        case 20:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QCustomPlot* >(); break;
            }
            break;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject PlotWindow::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_PlotWindow.data,
      qt_meta_data_PlotWindow,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *PlotWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *PlotWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_PlotWindow.stringdata0))
        return static_cast<void*>(this);
    return QWidget::qt_metacast(_clname);
}

int PlotWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 22)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 22;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 22)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 22;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
