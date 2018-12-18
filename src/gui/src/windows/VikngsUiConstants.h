#ifndef VIKNGSUICONSTANTS_H
#define VIKNGSUICONSTANTS_H
#include <QString>
#include <QColor>
#include <QPen>

const QString MAIN_LAYER = "main";
const QString OVERLAY_LAYER = "overlay";
const QString TOP_LAYER = "top";
const QString NO_LAYER = "dne";

const QColor BLACK = QColor(0,0,0, 255);
const QPen NO_PEN = Qt::NoPen;

inline QPen pen(QColor qcol, double alpha=1){
    return QPen(QColor(qcol.red(), qcol.green(), qcol.blue(), 255*alpha));
}

inline QBrush brush(QColor qcol, double alpha=1){
    return QBrush(QColor(qcol.red(), qcol.green(), qcol.blue(), 255*alpha));
}

#endif // VIKNGSUICONSTANTS_H
