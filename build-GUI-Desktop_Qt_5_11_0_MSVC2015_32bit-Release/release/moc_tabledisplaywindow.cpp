/****************************************************************************
** Meta object code from reading C++ file 'tabledisplaywindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.11.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../gui/src/windows/tabledisplaywindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'tabledisplaywindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.11.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_TableDisplayWindow_t {
    QByteArrayData data[15];
    char stringdata0[192];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_TableDisplayWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_TableDisplayWindow_t qt_meta_stringdata_TableDisplayWindow = {
    {
QT_MOC_LITERAL(0, 0, 18), // "TableDisplayWindow"
QT_MOC_LITERAL(1, 19, 10), // "initialize"
QT_MOC_LITERAL(2, 30, 0), // ""
QT_MOC_LITERAL(3, 31, 5), // "title"
QT_MOC_LITERAL(4, 37, 21), // "std::vector<Variant>&"
QT_MOC_LITERAL(5, 59, 8), // "variants"
QT_MOC_LITERAL(6, 68, 18), // "SimulationRequest&"
QT_MOC_LITERAL(7, 87, 7), // "request"
QT_MOC_LITERAL(8, 95, 5), // "index"
QT_MOC_LITERAL(9, 101, 17), // "fillGenotypeTable"
QT_MOC_LITERAL(10, 119, 12), // "variantIndex"
QT_MOC_LITERAL(11, 132, 16), // "fillVariantTable"
QT_MOC_LITERAL(12, 149, 31), // "on_table_variantTbl_cellClicked"
QT_MOC_LITERAL(13, 181, 3), // "row"
QT_MOC_LITERAL(14, 185, 6) // "column"

    },
    "TableDisplayWindow\0initialize\0\0title\0"
    "std::vector<Variant>&\0variants\0"
    "SimulationRequest&\0request\0index\0"
    "fillGenotypeTable\0variantIndex\0"
    "fillVariantTable\0on_table_variantTbl_cellClicked\0"
    "row\0column"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_TableDisplayWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    4,   34,    2, 0x0a /* Public */,
       9,    1,   43,    2, 0x08 /* Private */,
      11,    0,   46,    2, 0x08 /* Private */,
      12,    2,   47,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, QMetaType::QString, 0x80000000 | 4, 0x80000000 | 6, QMetaType::Int,    3,    5,    7,    8,
    QMetaType::Void, QMetaType::Int,   10,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,   13,   14,

       0        // eod
};

void TableDisplayWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        TableDisplayWindow *_t = static_cast<TableDisplayWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->initialize((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< std::vector<Variant>(*)>(_a[2])),(*reinterpret_cast< SimulationRequest(*)>(_a[3])),(*reinterpret_cast< int(*)>(_a[4]))); break;
        case 1: _t->fillGenotypeTable((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->fillVariantTable(); break;
        case 3: _t->on_table_variantTbl_cellClicked((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject TableDisplayWindow::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_TableDisplayWindow.data,
      qt_meta_data_TableDisplayWindow,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *TableDisplayWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *TableDisplayWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_TableDisplayWindow.stringdata0))
        return static_cast<void*>(this);
    return QWidget::qt_metacast(_clname);
}

int TableDisplayWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
