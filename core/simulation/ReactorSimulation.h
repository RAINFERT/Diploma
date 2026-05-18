#pragma once

#include "../dae/ReducedReactorDae.h"
#include "../models/WellStirredReactorModel.h"

#include <string>
#include <vector>

struct ReactorSimulationOptions
{
    double timeStepS = 1.0;
    int stepCount = 30;
};

struct ReactorSimulationPoint
{
    double timeS = 0.0;

    ReactorState state;
    WellStirredReactorEvaluation evaluation;
    MassBalanceResult massBalance;

    int newtonIterations = 0;
    double residualNorm = 0.0;

    bool converged = true;
    std::string message;
};

struct ReactorSimulationResult
{
    bool converged = true;
    std::string message;

    std::vector<ReactorSimulationPoint> points;

    const ReactorState& finalState() const;
};

class ReactorSimulationRunner
{
public:
    ReactorSimulationRunner(
        const ReducedReactorDae& dae,
        const WellStirredReactorModel& reactorModel
        );

    ReactorSimulationResult runRadauIIA2(
        const ReactorState& initialState,
        const ReactorSimulationOptions& options
        ) const;

    ReactorSimulationResult runRadauIIA3(
        const ReactorState& initialState,
        const ReactorSimulationOptions& options
        ) const;

    std::string resultTableToString(
        const ReactorSimulationResult& result,
        int printEvery = 1
        ) const;

private:
    const ReducedReactorDae& dae_;
    const WellStirredReactorModel& reactorModel_;

    ReactorSimulationPoint makePoint(
        double timeS,
        const ReactorState& state,
        int newtonIterations,
        double residualNorm,
        bool converged,
        const std::string& message
        ) const;
};
