#pragma once

#include "../dae/EnergyReactorDae.h"
#include "../models/WellStirredReactorModel.h"

#include <string>
#include <vector>

struct EnergyReactorSimulationOptions
{
    double timeStepS = 1.0;
    int stepCount = 30;
};

struct EnergyReactorSimulationPoint
{
    double timeS = 0.0;

    ReactorState state;

    WellStirredReactorEvaluation evaluation;
    MassBalanceResult massBalance;
    EnergyBalanceResult energyBalance;

    int newtonIterations = 0;
    double residualNorm = 0.0;

    bool converged = true;
    std::string message;
};

struct EnergyReactorSimulationResult
{
    bool converged = true;
    std::string message;

    std::vector<EnergyReactorSimulationPoint> points;

    const ReactorState& finalState() const;
};

class EnergyReactorSimulationRunner
{
public:
    EnergyReactorSimulationRunner(
        const EnergyReactorDae& dae,
        const WellStirredReactorModel& reactorModel
        );

    EnergyReactorSimulationResult runRadauIIA3(
        const ReactorState& initialState,
        const EnergyReactorSimulationOptions& options
        ) const;

    std::string resultTableToString(
        const EnergyReactorSimulationResult& result,
        int printEvery = 1
        ) const;

private:
    const EnergyReactorDae& dae_;
    const WellStirredReactorModel& reactorModel_;

    EnergyReactorSimulationPoint makePoint(
        double timeS,
        const ReactorState& state,
        int newtonIterations,
        double residualNorm,
        bool converged,
        const std::string& message
        ) const;
};
