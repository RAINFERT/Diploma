#pragma once

#include "../models/WellStirredReactorModel.h"
#include "../numerics/NewtonSolver.h"

#include <array>
#include <string>

struct ReducedReactorDaeVariables
{
    std::array<double, ComponentCount> massesKg{};

    double pressureBar = 0.0;
    double temperatureC = 0.0;
};

struct ReducedReactorDaeDerivatives
{
    std::array<double, ComponentCount> massDerivativesKgPerS{};
};

struct ReducedReactorDaeResidual
{
    std::array<double, ComponentCount> massResidualsKgPerS{};

    double volumeResidualM3 = 0.0;

    MassBalanceResult massBalance;
};

struct ReducedReactorDaeStepResult
{
    bool converged = false;

    int iterations = 0;

    double timeStepS = 0.0;

    ReactorState state;

    numerics::Vector variables;
    numerics::Vector residual;

    double residualNorm = 0.0;

    std::string method;
    std::string message;
};

class ReducedReactorDae
{
public:
    static constexpr std::size_t VariableCount =
        ComponentCount + 1;

    static constexpr std::size_t DerivativeCount =
        ComponentCount;

    static constexpr std::size_t ResidualCount =
        ComponentCount + 1;

    explicit ReducedReactorDae(
        const WellStirredReactorModel& reactorModel
        );

    ReducedReactorDaeStepResult implicitEulerStep(
        const ReactorState& oldState,
        double timeStepS
        ) const;

    ReducedReactorDaeStepResult radauIIA2Step(
        const ReactorState& oldState,
        double timeStepS
        ) const;

    ReducedReactorDaeStepResult radauIIA3Step(
        const ReactorState& oldState,
        double timeStepS
        ) const;

    std::string stepResultToString(
        const ReducedReactorDaeStepResult& result
        ) const;

    ReducedReactorDaeVariables variablesFromState(
        const ReactorState& state
        ) const;

    ReactorState stateFromVariables(
        const ReducedReactorDaeVariables& variables
        ) const;

    ReducedReactorDaeResidual computeResidual(
        const ReducedReactorDaeVariables& variables,
        const ReducedReactorDaeDerivatives& derivatives
        ) const;

    numerics::Vector variablesToVector(
        const ReducedReactorDaeVariables& variables
        ) const;

    ReducedReactorDaeVariables variablesFromVector(
        const numerics::Vector& vector,
        double temperatureC
        ) const;

    numerics::Vector derivativesToVector(
        const ReducedReactorDaeDerivatives& derivatives
        ) const;

    ReducedReactorDaeDerivatives derivativesFromVector(
        const numerics::Vector& vector
        ) const;

    numerics::Vector residualToVector(
        const ReducedReactorDaeResidual& residual
        ) const;

    numerics::Vector computeResidualVector(
        const numerics::Vector& variablesVector,
        const numerics::Vector& derivativesVector,
        double temperatureC
        ) const;

    std::string residualToString(
        const ReducedReactorDaeResidual& residual
        ) const;

private:
    const WellStirredReactorModel& reactorModel_;
};
