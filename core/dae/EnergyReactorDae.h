#pragma once

#include "../models/WellStirredReactorModel.h"
#include "../numerics/NewtonSolver.h"

#include <array>
#include <cstddef>
#include <string>

struct EnergyReactorDaeVariables
{
    std::array<double, ComponentCount> massesKg{};

    double energyJ = 0.0;

    double pressureBar = 0.0;
    double temperatureC = 0.0;
};

struct EnergyReactorDaeDerivatives
{
    std::array<double, ComponentCount> massDerivativesKgPerS{};

    double energyDerivativeW = 0.0;
};

struct EnergyReactorDaeResidual
{
    std::array<double, ComponentCount> massResidualsKgPerS{};

    double energyResidualW = 0.0;

    double volumeResidualM3 = 0.0;

    double inventoryEnergyResidualJ = 0.0;

    MassBalanceResult massBalance;
    EnergyBalanceResult energyBalance;
};

struct EnergyReactorDaeStepResult
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

class EnergyReactorDae
{
public:
    static constexpr std::size_t VariableCount =
        ComponentCount + 3;
    // masses + E + P + T

    static constexpr std::size_t DerivativeCount =
        ComponentCount + 1;
    // dM/dt + dE/dt

    static constexpr std::size_t ResidualCount =
        ComponentCount + 3;
    // R_M + R_E + R_V + R_H

    explicit EnergyReactorDae(
        const WellStirredReactorModel& reactorModel
        );

    EnergyReactorDaeVariables variablesFromState(
        const ReactorState& state
        ) const;

    ReactorState stateFromVariables(
        const EnergyReactorDaeVariables& variables
        ) const;

    EnergyReactorDaeResidual computeResidual(
        const EnergyReactorDaeVariables& variables,
        const EnergyReactorDaeDerivatives& derivatives
        ) const;

    numerics::Vector variablesToVector(
        const EnergyReactorDaeVariables& variables
        ) const;

    EnergyReactorDaeVariables variablesFromVector(
        const numerics::Vector& vector
        ) const;

    numerics::Vector derivativesToVector(
        const EnergyReactorDaeDerivatives& derivatives
        ) const;

    EnergyReactorDaeDerivatives derivativesFromVector(
        const numerics::Vector& vector
        ) const;

    numerics::Vector residualToVector(
        const EnergyReactorDaeResidual& residual
        ) const;

    numerics::Vector computeResidualVector(
        const numerics::Vector& variablesVector,
        const numerics::Vector& derivativesVector
        ) const;

    std::string residualToString(
        const EnergyReactorDaeResidual& residual
        ) const;

    EnergyReactorDaeStepResult radauIIA3Step(
        const ReactorState& oldState,
        double timeStepS
        ) const;

    std::string stepResultToString(
        const EnergyReactorDaeStepResult& result
        ) const;

private:
    const WellStirredReactorModel& reactorModel_;
};
