class Runner : public QObject {
    Q_OBJECT

public:
    Runner();
    ~Runner();

public slots:
    void process();

signals:
    void finished();

private:
    // add your variables here
};


