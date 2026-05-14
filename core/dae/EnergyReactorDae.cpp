#include "EnergyReactorDae.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>

EnergyReactorDae::EnergyReactorDae(
    const WellStirredReactorModel& reactorModel
    )
    : reactorModel_(reactorModel)
{
}

EnergyReactorDaeVariables EnergyReactorDae::variablesFromState(
    const ReactorState& state
    ) const
{
    EnergyReactorDaeVariables variables;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const Component component =
            static_cast<Component>(i);

        variables.massesKg[i] =
            state.massKg(component);
    }

    variables.energyJ =
        state.energyJ();

    variables.pressureBar =
        state.pressureBar();

    variables.temperatureC =
        state.temperatureC();

    return variables;
}

ReactorState EnergyReactorDae::stateFromVariables(
    const EnergyReactorDaeVariables& variables
    ) const
{
    ReactorState state;

    double massScaleKg = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        massScaleKg +=
            std::abs(variables.massesKg[i]);
    }

    const double negativeMassToleranceKg =
        1.0e-12
        + 1.0e-8 * std::max(1.0, massScaleKg);

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        double massKg =
            variables.massesKg[i];

        if (massKg < -negativeMassToleranceKg)
        {
            throw std::runtime_error("Energy DAE variable mass cannot be negative");
        }

        if (massKg < 0.0)
        {
            massKg = 0.0;
        }

        const Component component =
            static_cast<Component>(i);

        state.setMassKg(
            component,
            massKg
            );
    }

    if (variables.pressureBar <= 0.0)
    {
        throw std::runtime_error("Energy DAE pressure must be positive");
    }

    if (variables.temperatureC <= -273.15)
    {
        throw std::runtime_error("Energy DAE temperature is below absolute zero");
    }

    state.setEnergyJ(
        variables.energyJ
        );

    state.setPressureBar(
        variables.pressureBar
        );

    state.setTemperatureC(
        variables.temperatureC
        );

    return state;
}

EnergyReactorDaeResidual EnergyReactorDae::computeResidual(
    const EnergyReactorDaeVariables& variables,
    const EnergyReactorDaeDerivatives& derivatives
    ) const
{
    EnergyReactorDaeResidual residual;

    const ReactorState state =
        stateFromVariables(
            variables
            );

    residual.massBalance =
        reactorModel_.computeMassBalance(
            state
            );

    residual.energyBalance =
        reactorModel_.computeEnergyBalance(
            state
            );

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        residual.massResidualsKgPerS[i] =
            derivatives.massDerivativesKgPerS[i]
            - residual.massBalance.massDerivativesKgPerS[i];
    }

    residual.energyResidualW =
        derivatives.energyDerivativeW
        - residual.energyBalance.energyDerivativeW;

    residual.volumeResidualM3 =
        residual.massBalance.evaluation.volumeResidualM3;

    residual.inventoryEnergyResidualJ =
        residual.energyBalance.energyResidualJ;

    return residual;
}

numerics::Vector EnergyReactorDae::variablesToVector(
    const EnergyReactorDaeVariables& variables
    ) const
{
    numerics::Vector vector(
        VariableCount,
        0.0
        );

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        vector[i] =
            variables.massesKg[i];
    }

    vector[ComponentCount] =
        variables.energyJ;

    vector[ComponentCount + 1] =
        variables.pressureBar;

    vector[ComponentCount + 2] =
        variables.temperatureC;

    return vector;
}

EnergyReactorDaeVariables EnergyReactorDae::variablesFromVector(
    const numerics::Vector& vector
    ) const
{
    if (vector.size() != VariableCount)
    {
        throw std::invalid_argument("Energy DAE variable vector has invalid size");
    }

    EnergyReactorDaeVariables variables;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        variables.massesKg[i] =
            vector[i];
    }

    variables.energyJ =
        vector[ComponentCount];

    variables.pressureBar =
        vector[ComponentCount + 1];

    variables.temperatureC =
        vector[ComponentCount + 2];

    return variables;
}

numerics::Vector EnergyReactorDae::derivativesToVector(
    const EnergyReactorDaeDerivatives& derivatives
    ) const
{
    numerics::Vector vector(
        DerivativeCount,
        0.0
        );

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        vector[i] =
            derivatives.massDerivativesKgPerS[i];
    }

    vector[ComponentCount] =
        derivatives.energyDerivativeW;

    return vector;
}

EnergyReactorDaeDerivatives EnergyReactorDae::derivativesFromVector(
    const numerics::Vector& vector
    ) const
{
    if (vector.size() != DerivativeCount)
    {
        throw std::invalid_argument("Energy DAE derivative vector has invalid size");
    }

    EnergyReactorDaeDerivatives derivatives;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        derivatives.massDerivativesKgPerS[i] =
            vector[i];
    }

    derivatives.energyDerivativeW =
        vector[ComponentCount];

    return derivatives;
}

numerics::Vector EnergyReactorDae::residualToVector(
    const EnergyReactorDaeResidual& residual
    ) const
{
    numerics::Vector vector(
        ResidualCount,
        0.0
        );

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        vector[i] =
            residual.massResidualsKgPerS[i];
    }

    vector[ComponentCount] =
        residual.energyResidualW;

    vector[ComponentCount + 1] =
        residual.volumeResidualM3;

    vector[ComponentCount + 2] =
        residual.inventoryEnergyResidualJ;

    return vector;
}

numerics::Vector EnergyReactorDae::computeResidualVector(
    const numerics::Vector& variablesVector,
    const numerics::Vector& derivativesVector
    ) const
{
    const EnergyReactorDaeVariables variables =
        variablesFromVector(
            variablesVector
            );

    const EnergyReactorDaeDerivatives derivatives =
        derivativesFromVector(
            derivativesVector
            );

    const EnergyReactorDaeResidual residual =
        computeResidual(
            variables,
            derivatives
            );

    return residualToVector(
        residual
        );
}

std::string EnergyReactorDae::residualToString(
    const EnergyReactorDaeResidual& residual
    ) const
{
    std::ostringstream out;

    out << std::fixed << std::setprecision(10);

    out << "Energy reactor DAE residual:\n";

    out << "  Mass residuals:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "    R_M[" << i << "] = "
            << residual.massResidualsKgPerS[i]
            << " kg/s\n";
    }

    out << "  Energy residual:\n";
    out << "    R_E = "
        << residual.energyResidualW
        << " W\n";

    out << "  Algebraic residuals:\n";
    out << "    R_V = "
        << residual.volumeResidualM3
        << " m3\n";
    out << "    R_H = "
        << residual.inventoryEnergyResidualJ
        << " J\n";

    out << "  Underlying balances:\n";

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "    f_M[" << i << "] = "
            << residual.massBalance.massDerivativesKgPerS[i]
            << " kg/s\n";
    }

    out << "    f_E = "
        << residual.energyBalance.energyDerivativeW
        << " W\n";

    out << "  State diagnostics:\n";
    out << "    P      = "
        << residual.massBalance.evaluation.pressurePa / 1.0e5
        << " bar\n";
    out << "    T      = "
        << residual.massBalance.evaluation.temperatureK - 273.15
        << " C\n";
    out << "    beta   = "
        << residual.massBalance.evaluation.flash.beta
        << "\n";
    out << "    V_calc = "
        << residual.massBalance.evaluation.calculatedVolumeM3
        << " m3\n";
    out << "    E_inv  = "
        << residual.energyBalance.inventoryEnergyJ
        << " J\n";

    return out.str();
}

EnergyReactorDaeStepResult EnergyReactorDae::radauIIA3Step(
    const ReactorState& oldState,
    double timeStepS
    ) const
{
    if (timeStepS <= 0.0)
    {
        throw std::invalid_argument("Time step must be positive");
    }

    constexpr int StageCount = 3;

    // Stage unknowns:
    // dM(C2H6)/dt, dM(C5H12)/dt, dM(H2O)/dt, dE/dt, P, T
    constexpr std::size_t StageUnknownCount =
        ComponentCount + 3;

    const double sqrt6 =
        std::sqrt(6.0);

    const double c1 =
        (4.0 - sqrt6) / 10.0;

    const double c2 =
        (4.0 + sqrt6) / 10.0;

    const double c3 =
        1.0;

    const double a11 =
        (88.0 - 7.0 * sqrt6) / 360.0;

    const double a12 =
        (296.0 - 169.0 * sqrt6) / 1800.0;

    const double a13 =
        (-2.0 + 3.0 * sqrt6) / 225.0;

    const double a21 =
        (296.0 + 169.0 * sqrt6) / 1800.0;

    const double a22 =
        (88.0 + 7.0 * sqrt6) / 360.0;

    const double a23 =
        (-2.0 - 3.0 * sqrt6) / 225.0;

    const double a31 =
        (16.0 - sqrt6) / 36.0;

    const double a32 =
        (16.0 + sqrt6) / 36.0;

    const double a33 =
        1.0 / 9.0;

    const EnergyReactorDaeVariables oldVariables =
        variablesFromState(oldState);

    const numerics::Vector oldVariableVector =
        variablesToVector(oldVariables);

    const MassBalanceResult oldMassBalance =
        reactorModel_.computeMassBalance(oldState);

    const EnergyBalanceResult oldEnergyBalance =
        reactorModel_.computeEnergyBalance(oldState);

    auto makePressureGuess =
        [this, &oldState](const std::array<double, ComponentCount>& masses)
    {
        ReactorState guessState =
            oldState;

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            const Component component =
                static_cast<Component>(i);

            guessState.setMassKg(
                component,
                std::max(0.0, masses[i])
                );
        }

        const PressureInitializationResult pressureGuess =
            reactorModel_.initializePressureForVolumeNewton(
                guessState,
                1.0e4,
                1.0e7,
                1.0e-6,
                15,
                true
                );

        if (pressureGuess.converged)
        {
            return pressureGuess.pressureBar;
        }

        return oldState.pressureBar();
    };

    numerics::Vector initialGuess(
        StageCount * StageUnknownCount,
        0.0
        );

    for (int stage = 0; stage < StageCount; ++stage)
    {
        double c = c3;

        if (stage == 0)
        {
            c = c1;
        }
        else if (stage == 1)
        {
            c = c2;
        }

        std::array<double, ComponentCount> stageMassGuess{};

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            stageMassGuess[i] =
                oldVariables.massesKg[i]
                + timeStepS
                      * c
                      * oldMassBalance.massDerivativesKgPerS[i];

            initialGuess[
                static_cast<std::size_t>(stage) * StageUnknownCount
                + i
            ] =
                oldMassBalance.massDerivativesKgPerS[i];
        }

        initialGuess[
            static_cast<std::size_t>(stage) * StageUnknownCount
            + ComponentCount
        ] =
            oldEnergyBalance.energyDerivativeW;

        initialGuess[
            static_cast<std::size_t>(stage) * StageUnknownCount
            + ComponentCount + 1
        ] =
            makePressureGuess(stageMassGuess);

        initialGuess[
            static_cast<std::size_t>(stage) * StageUnknownCount
            + ComponentCount + 2
        ] =
            oldState.temperatureC();
    }

    auto stageMassDerivative =
        [](const numerics::Vector& unknowns, int stage, std::size_t component)
    {
        constexpr std::size_t blockSize =
            ComponentCount + 3;

        return unknowns[
            static_cast<std::size_t>(stage) * blockSize
            + component
        ];
    };

    auto stageEnergyDerivative =
        [](const numerics::Vector& unknowns, int stage)
    {
        constexpr std::size_t blockSize =
            ComponentCount + 3;

        return unknowns[
            static_cast<std::size_t>(stage) * blockSize
            + ComponentCount
        ];
    };

    auto stagePressure =
        [](const numerics::Vector& unknowns, int stage)
    {
        constexpr std::size_t blockSize =
            ComponentCount + 3;

        return unknowns[
            static_cast<std::size_t>(stage) * blockSize
            + ComponentCount + 1
        ];
    };

    auto stageTemperature =
        [](const numerics::Vector& unknowns, int stage)
    {
        constexpr std::size_t blockSize =
            ComponentCount + 3;

        return unknowns[
            static_cast<std::size_t>(stage) * blockSize
            + ComponentCount + 2
        ];
    };

    double massResidualScaleKgPerS = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        massResidualScaleKgPerS +=
            std::abs(oldMassBalance.massDerivativesKgPerS[i]);
    }

    massResidualScaleKgPerS =
        std::max(
            1.0e-3,
            massResidualScaleKgPerS
            );

    const double energyResidualScaleW =
        std::max(
            1.0,
            std::abs(oldEnergyBalance.inletEnthalpyFlowW)
                + std::abs(oldEnergyBalance.outletEnthalpyFlowW)
                + std::abs(oldEnergyBalance.heatTransferW)
                + std::abs(oldEnergyBalance.energyDerivativeW)
            );

    const double volumeResidualScaleM3 =
        std::max(
            1.0,
            reactorModel_.parameters().volumeM3
            );

    const double inventoryEnergyResidualScaleJ =
        std::max(
            1.0e6,
            std::abs(oldVariables.energyJ)
            );

    auto residualFunction =
        [this,
         oldVariables,
         timeStepS,
         stageMassDerivative,
         stageEnergyDerivative,
         stagePressure,
         stageTemperature,
         massResidualScaleKgPerS,
         energyResidualScaleW,
         volumeResidualScaleM3,
         inventoryEnergyResidualScaleJ,
         a11,
         a12,
         a13,
         a21,
         a22,
         a23,
         a31,
         a32,
         a33]
        (const numerics::Vector& unknowns)
    {
        constexpr std::size_t totalUnknownCount =
            3 * (ComponentCount + 3);

        if (unknowns.size() != totalUnknownCount)
        {
            throw std::runtime_error("Energy Radau IIA 3 unknown vector has invalid size");
        }

        numerics::Vector totalResidual(
            totalUnknownCount,
            0.0
            );

        for (int stage = 0; stage < 3; ++stage)
        {
            EnergyReactorDaeVariables variables;
            EnergyReactorDaeDerivatives derivatives;

            for (std::size_t i = 0; i < ComponentCount; ++i)
            {
                const double k1 =
                    stageMassDerivative(unknowns, 0, i);

                const double k2 =
                    stageMassDerivative(unknowns, 1, i);

                const double k3 =
                    stageMassDerivative(unknowns, 2, i);

                double mass =
                    oldVariables.massesKg[i];

                if (stage == 0)
                {
                    mass +=
                        timeStepS
                        * (
                            a11 * k1
                            + a12 * k2
                            + a13 * k3
                            );
                }
                else if (stage == 1)
                {
                    mass +=
                        timeStepS
                        * (
                            a21 * k1
                            + a22 * k2
                            + a23 * k3
                            );
                }
                else
                {
                    mass +=
                        timeStepS
                        * (
                            a31 * k1
                            + a32 * k2
                            + a33 * k3
                            );
                }

                variables.massesKg[i] =
                    mass;

                derivatives.massDerivativesKgPerS[i] =
                    stageMassDerivative(
                        unknowns,
                        stage,
                        i
                        );
            }

            const double e1 =
                stageEnergyDerivative(unknowns, 0);

            const double e2 =
                stageEnergyDerivative(unknowns, 1);

            const double e3 =
                stageEnergyDerivative(unknowns, 2);

            double energy =
                oldVariables.energyJ;

            if (stage == 0)
            {
                energy +=
                    timeStepS
                    * (
                        a11 * e1
                        + a12 * e2
                        + a13 * e3
                        );
            }
            else if (stage == 1)
            {
                energy +=
                    timeStepS
                    * (
                        a21 * e1
                        + a22 * e2
                        + a23 * e3
                        );
            }
            else
            {
                energy +=
                    timeStepS
                    * (
                        a31 * e1
                        + a32 * e2
                        + a33 * e3
                        );
            }

            variables.energyJ =
                energy;

            variables.pressureBar =
                stagePressure(
                    unknowns,
                    stage
                    );

            variables.temperatureC =
                stageTemperature(
                    unknowns,
                    stage
                    );

            derivatives.energyDerivativeW =
                stageEnergyDerivative(
                    unknowns,
                    stage
                    );

            const EnergyReactorDaeResidual stageResidual =
                computeResidual(
                    variables,
                    derivatives
                    );

            const numerics::Vector stageResidualVector =
                residualToVector(
                    stageResidual
                    );

            const std::size_t offset =
                static_cast<std::size_t>(stage) * ResidualCount;

            for (std::size_t i = 0; i < ComponentCount; ++i)
            {
                totalResidual[offset + i] =
                    stageResidualVector[i]
                    / massResidualScaleKgPerS;
            }

            totalResidual[offset + ComponentCount] =
                stageResidualVector[ComponentCount]
                / energyResidualScaleW;

            totalResidual[offset + ComponentCount + 1] =
                stageResidualVector[ComponentCount + 1]
                / volumeResidualScaleM3;

            totalResidual[offset + ComponentCount + 2] =
                stageResidualVector[ComponentCount + 2]
                / inventoryEnergyResidualScaleJ;
        }

        return totalResidual;
    };

    numerics::NewtonOptions options;

    options.maxIterations = 80;
    options.residualTolerance = 1.0e-8;
    options.stepTolerance = 1.0e-10;
    options.finiteDifferenceStep = 1.0e-5;
    options.minDampingFactor = 1.0e-8;

    const double inf =
        std::numeric_limits<double>::infinity();

    options.lowerBounds =
        numerics::Vector{
            -inf, -inf, -inf, -inf, 0.01, -273.14,
            -inf, -inf, -inf, -inf, 0.01, -273.14,
            -inf, -inf, -inf, -inf, 0.01, -273.14
        };

    const numerics::NewtonResult newtonResult =
        numerics::NewtonSolver::solve(
            residualFunction,
            initialGuess,
            options
            );

    EnergyReactorDaeStepResult result;

    result.converged =
        newtonResult.converged;

    result.iterations =
        newtonResult.iterations;

    result.timeStepS =
        timeStepS;

    result.residual =
        newtonResult.residual;

    result.residualNorm =
        newtonResult.residualNorm;

    result.message =
        newtonResult.message;

    result.method =
        "energy_radau_iia_3";

    EnergyReactorDaeVariables finalVariables;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const double k1 =
            stageMassDerivative(
                newtonResult.solution,
                0,
                i
                );

        const double k2 =
            stageMassDerivative(
                newtonResult.solution,
                1,
                i
                );

        const double k3 =
            stageMassDerivative(
                newtonResult.solution,
                2,
                i
                );

        finalVariables.massesKg[i] =
            oldVariables.massesKg[i]
            + timeStepS
                  * (
                      a31 * k1
                      + a32 * k2
                      + a33 * k3
                      );
    }

    const double e1 =
        stageEnergyDerivative(
            newtonResult.solution,
            0
            );

    const double e2 =
        stageEnergyDerivative(
            newtonResult.solution,
            1
            );

    const double e3 =
        stageEnergyDerivative(
            newtonResult.solution,
            2
            );

    finalVariables.energyJ =
        oldVariables.energyJ
        + timeStepS
              * (
                  a31 * e1
                  + a32 * e2
                  + a33 * e3
                  );

    finalVariables.pressureBar =
        stagePressure(
            newtonResult.solution,
            2
            );

    finalVariables.temperatureC =
        stageTemperature(
            newtonResult.solution,
            2
            );

    result.state =
        stateFromVariables(
            finalVariables
            );

    result.variables =
        variablesToVector(
            finalVariables
            );

    return result;
}

std::string EnergyReactorDae::stepResultToString(
    const EnergyReactorDaeStepResult& result
    ) const
{
    std::ostringstream out;

    out << std::fixed << std::setprecision(10);

    out << "Energy DAE step result:\n";
    out << "  method        = "
        << result.method << "\n";
    out << "  converged     = "
        << (result.converged ? "true" : "false") << "\n";
    out << "  iterations    = "
        << result.iterations << "\n";
    out << "  dt            = "
        << result.timeStepS << " s\n";
    out << "  residual norm = "
        << result.residualNorm << "\n";
    out << "  message       = "
        << result.message << "\n";

    out << "\n  State after step:\n";
    out << result.state.toString();

    out << "\n  Residual vector:\n";

    for (std::size_t i = 0; i < result.residual.size(); ++i)
    {
        out << "    R[" << i << "] = "
            << result.residual[i] << "\n";
    }

    return out.str();
}
