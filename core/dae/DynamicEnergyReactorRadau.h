#pragma once

#include "DynamicEnergyReactorDae.h"
#include "../models/DynamicWellStirredReactorModel.h"
#include "../numerics/NewtonSolver.h"

#include <string>
#include <vector>

struct DynamicRadauOptions {
    double residualTolerance = 1.0e-8;
    double stepTolerance = 1.0e-10;
    int maxNewtonIterations = 30;
};

struct DynamicRadauStepResult {
    bool converged = false;
    int iterations = 0;
    double residualNorm = 0.0;

    std::string message;

    DynamicVector oldVariables;
    DynamicVector newVariables;
};

class DynamicEnergyReactorRadau {
public:
    DynamicEnergyReactorRadau(
        const DynamicEnergyReactorDae& dae,
        const DynamicWellStirredReactorModel& reactorModel
        );

    DynamicRadauStepResult radauIIA3Step(
        const DynamicVector& oldVariables,
        double timeStepS,
        const DynamicRadauOptions& options = DynamicRadauOptions{}
        ) const;

private:
    const DynamicEnergyReactorDae& dae_;
    const DynamicWellStirredReactorModel& reactorModel_;

    std::size_t differentialCount() const;

    std::size_t differentialVariableIndex(
        std::size_t differentialIndex
        ) const;

    std::size_t stageBlockSize() const;

    std::size_t stageOffset(
        std::size_t stageIndex
        ) const;

    std::size_t stageVariableIndex(
        std::size_t stageIndex,
        std::size_t variableIndex
        ) const;

    std::size_t stageDerivativeIndex(
        std::size_t stageIndex,
        std::size_t differentialIndex
        ) const;

    DynamicVector buildInitialGuess(
        const DynamicVector& oldVariables,
        double timeStepS
        ) const;

    DynamicVector extractStageVariables(
        const DynamicVector& unknowns,
        std::size_t stageIndex
        ) const;

    DynamicVector extractStageDerivativesFull(
        const DynamicVector& unknowns,
        std::size_t stageIndex
        ) const;

    DynamicVector buildResidual(
        const DynamicVector& unknowns,
        const DynamicVector& oldVariables,
        double timeStepS
        ) const;

    DynamicVector buildNewVariables(
        const DynamicVector& oldVariables,
        const DynamicVector& solvedUnknowns,
        double timeStepS
        ) const;
};
