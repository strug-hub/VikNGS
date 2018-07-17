/****************************************************************************
** Meta object code from reading C++ file 'runner.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.11.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../gui/src/runner.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#include <QtCore/QVector>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'runner.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.11.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_Runner_t {
    QByteArrayData data[14];
    char stringdata0[190];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Runner_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Runner_t qt_meta_stringdata_Runner = {
    {
QT_MOC_LITERAL(0, 0, 6), // "Runner"
QT_MOC_LITERAL(1, 7, 11), // "jobFinished"
QT_MOC_LITERAL(2, 19, 0), // ""
QT_MOC_LITERAL(3, 20, 16), // "QVector<Variant>"
QT_MOC_LITERAL(4, 37, 18), // "simulationFinished"
QT_MOC_LITERAL(5, 56, 34), // "std::vector<std::vector<Varia..."
QT_MOC_LITERAL(6, 91, 17), // "SimulationRequest"
QT_MOC_LITERAL(7, 109, 8), // "complete"
QT_MOC_LITERAL(8, 118, 9), // "runVikngs"
QT_MOC_LITERAL(9, 128, 13), // "runSimulation"
QT_MOC_LITERAL(10, 142, 10), // "setRequest"
QT_MOC_LITERAL(11, 153, 7), // "Request"
QT_MOC_LITERAL(12, 161, 7), // "request"
QT_MOC_LITERAL(13, 169, 20) // "setSimulationRequest"

    },
    "Runner\0jobFinished\0\0QVector<Variant>\0"
    "simulationFinished\0"
    "std::vector<std::vector<Variant> >\0"
    "SimulationRequest\0complete\0runVikngs\0"
    "runSimulation\0setRequest\0Request\0"
    "request\0setSimulationRequest"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Runner[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   49,    2, 0x06 /* Public */,
       4,    2,   52,    2, 0x06 /* Public */,
       7,    0,   57,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       8,    0,   58,    2, 0x0a /* Public */,
       9,    0,   59,    2, 0x0a /* Public */,
      10,    1,   60,    2, 0x0a /* Public */,
      13,    1,   63,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    2,
    QMetaType::Void, 0x80000000 | 5, 0x80000000 | 6,    2,    2,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 11,   12,
    QMetaType::Void, 0x80000000 | 6,   12,

       0        // eod
};

void Runner::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Runner *_t = static_cast<Runner *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->jobFinished((*reinterpret_cast< QVector<Variant>(*)>(_a[1]))); break;
        case 1: _t->simulationFinished((*reinterpret_cast< std::vector<std::vector<Variant> >(*)>(_a[1])),(*reinterpret_cast< SimulationRequest(*)>(_a[2]))); break;
        case 2: _t->complete(); break;
        case 3: _t->runVikngs(); break;
        case 4: _t->runSimulation(); break;
        case 5: _t->setRequest((*reinterpret_cast< Request(*)>(_a[1]))); break;
        case 6: _t->setSimulationRequest((*reinterpret_cast< SimulationRequest(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (Runner::*)(QVector<Variant> );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Runner::jobFinished)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (Runner::*)(std::vector<std::vector<Variant>> , SimulationRequest );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Runner::simulationFinished)) {
                *result = 1;
                return;
            }
        }
        {
            using _t = void (Runner::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Runner::complete)) {
                *result = 2;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject Runner::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Runner.data,
      qt_meta_data_Runner,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *Runner::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Runner::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_Runner.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int Runner::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}

// SIGNAL 0
void Runner::jobFinished(QVector<Variant> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void Runner::simulationFinished(std::vector<std::vector<Variant>> _t1, SimulationRequest _t2)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void Runner::complete()
{
    QMetaObject::activate(this, &staticMetaObject, 2, nullptr);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
