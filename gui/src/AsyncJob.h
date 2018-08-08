#include <QObject>
#include "../src/vikNGS.h"
#include "../src/Variant.h"
#include "simulation/simulation.h"

class AsyncJob : public QObject {
    Q_OBJECT

public:
    AsyncJob(){}
    ~AsyncJob(){}

public slots:
    void runVikngs() {
        try{
            Data result = startVikNGS(request);
            emit jobFinished(result);
        }
        catch(...){
            Data empty;
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
    void jobFinished(Data);
    void simulationFinished(std::vector<std::vector<Variant>>, SimulationRequest);
    void complete();

private:
    Request request;
    SimulationRequest simRequest;
};

