#include <QObject>
#include <QString>
#include <QColor>

class QLog : public QObject {
    Q_OBJECT

public:
    QLog(){}
    ~QLog(){}

public slots:
    void out(QString str, QColor c) {
        emit pushOutput(str, c);
    }

signals:
    void pushOutput(QString str, QColor c);

};

QLog* getQLog();
