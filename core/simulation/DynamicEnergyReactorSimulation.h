#pragma once

#include "../dae/DynamicEnergyReactorRadau.h"
#include "../dae/DynamicEnergyReactorDae.h"
#include "../models/DynamicWellStirredReactorModel.h"

#include <string>
#include <vector>

struct DynamicEnergyReactorSimulationOptions {
    double timeStepS = 0.1;
    int stepCount = 10;

    DynamicRadauOptions radauOptions;
};

struct DynamicEnergyReactorSimulationPoint {
    double timeS = 0.0;

    DynamicVector variables;
    DynamicEnergyReactorDaeEvaluation evaluation;

    DynamicMassBalanceResult massBalance;
    DynamicEnergyBalanceResult energyBalance;

    bool converged = true;
    int newtonIterations = 0;
    double residualNorm = 0.0;
};

struct DynamicEnergyReactorSimulationResult {
    bool converged = true;
    std::string message;

    std::vector<DynamicEnergyReactorSimulationPoint> points;

    const DynamicEnergyReactorSimulationPoint& finalPoint() const {
        return points.back();
    }
};

class DynamicEnergyReactorSimulationRunner {
public:
    DynamicEnergyReactorSimulationRunner(
        const DynamicEnergyReactorDae& dae,
        const DynamicWellStirredReactorModel& reactorModel
        );

    DynamicEnergyReactorSimulationResult runRadauIIA3(
        const DynamicVector& initialVariables,
        const DynamicEnergyReactorSimulationOptions& options
        ) const;

    std::string resultTableToString(
        const DynamicEnergyReactorSimulationResult& result,
        const ComponentSet& components,
        int printEvery = 1
        ) const;

private:
    const DynamicEnergyReactorDae& dae_;
    const DynamicWellStirredReactorModel& reactorModel_;
};
