#include <QObject>
#include <QVector>
#include "../src/RVS.h"
#include "../src/Variant.h"
#include "simulation/simulation.h"

class Runner : public QObject {
    Q_OBJECT

public:
    Runner(){}
    ~Runner(){}

public slots:
    void runVikngs() {
        try{
            std::vector<Variant> v = startVikNGS(request);
            QVector<Variant> variants = QVector<Variant>::fromStdVector(v);
            emit jobFinished(variants);
        }
        catch(...){
            QVector<Variant> empty;
            emit jobFinished(empty);
        }

        emit complete();
    }

    void runSimulation() {
        try{
            std::vector<std::vector<Variant>> results = startSimulation(simRequest);
            emit simulationFinished(results, simRequest);
        }
        catch(...){
            std::vector<std::vector<Variant>> empty;
            emit simulationFinished(empty, simRequest);
        }

        emit complete();
    }

    void setRequest(Request request){
        this->request = request;
    }
    void setSimulationRequest(SimulationRequest request){
        this->simRequest = request;
    }

signals:
    void jobFinished(QVector<Variant>);
    void simulationFinished(std::vector<std::vector<Variant>>, SimulationRequest);
    void complete();

private:
    Request request;
    SimulationRequest simRequest;
};

