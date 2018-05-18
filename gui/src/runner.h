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
            std::vector<std::vector<Variant>> results = startSimulation(simRequests);
            emit simulationFinished(results, simRequests);
        }
        catch(...){
            std::vector<std::vector<Variant>> empty;
            emit simulationFinished(empty, simRequests);
        }

        emit complete();
    }

    void setRequest(Request request){
        this->request = request;
    }
    void setSimulationRequests( std::vector<SimulationRequest> requests){
        this->simRequests = requests;
    }

signals:
    void jobFinished(QVector<Variant>);
    void simulationFinished(std::vector<std::vector<Variant>>, std::vector<SimulationRequest>);
    void complete();

private:
    Request request;
    std::vector<SimulationRequest> simRequests;
};

