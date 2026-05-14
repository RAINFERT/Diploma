#include "ReducedReactorDae.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <limits>

ReducedReactorDae::ReducedReactorDae(
    const WellStirredReactorModel& reactorModel
    )
    : reactorModel_(reactorModel)
{
}

ReducedReactorDaeVariables ReducedReactorDae::variablesFromState(
    const ReactorState& state
    ) const
{
    ReducedReactorDaeVariables variables;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const Component component =
            static_cast<Component>(i);

        variables.massesKg[i] =
            state.massKg(component);
    }

    variables.pressureBar =
        state.pressureBar();

    variables.temperatureC =
        state.temperatureC();

    return variables;
}

ReactorState ReducedReactorDae::stateFromVariables(
    const ReducedReactorDaeVariables& variables
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
            throw std::runtime_error("DAE variable mass cannot be negative");
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
        throw std::runtime_error("DAE pressure must be positive");
    }

    if (variables.temperatureC <= -273.15)
    {
        throw std::runtime_error("DAE temperature is below absolute zero");
    }

    state.setPressureBar(
        variables.pressureBar
        );

    state.setTemperatureC(
        variables.temperatureC
        );

    return state;
}

ReducedReactorDaeResidual ReducedReactorDae::computeResidual(
    const ReducedReactorDaeVariables& variables,
    const ReducedReactorDaeDerivatives& derivatives
    ) const
{
    ReducedReactorDaeResidual residual;

    ReactorState state =
        stateFromVariables(variables);

    residual.massBalance =
        reactorModel_.computeMassBalance(state);

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        residual.massResidualsKgPerS[i] =
            derivatives.massDerivativesKgPerS[i]
            - residual.massBalance.massDerivativesKgPerS[i];
    }

    residual.volumeResidualM3 =
        residual.massBalance.evaluation.volumeResidualM3;

    return residual;
}

numerics::Vector ReducedReactorDae::variablesToVector(
    const ReducedReactorDaeVariables& variables
    ) const
{
    numerics::Vector vector(VariableCount, 0.0);

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        vector[i] = variables.massesKg[i];
    }

    vector[ComponentCount] =
        variables.pressureBar;

    return vector;
}

ReducedReactorDaeVariables ReducedReactorDae::variablesFromVector(
    const numerics::Vector& vector,
    double temperatureC
    ) const
{
    if (vector.size() != VariableCount)
    {
        throw std::invalid_argument("Reduced DAE variable vector has invalid size");
    }

    ReducedReactorDaeVariables variables;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        variables.massesKg[i] = vector[i];
    }

    variables.pressureBar =
        vector[ComponentCount];

    variables.temperatureC =
        temperatureC;

    return variables;
}

numerics::Vector ReducedReactorDae::derivativesToVector(
    const ReducedReactorDaeDerivatives& derivatives
    ) const
{
    numerics::Vector vector(DerivativeCount, 0.0);

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        vector[i] =
            derivatives.massDerivativesKgPerS[i];
    }

    return vector;
}

ReducedReactorDaeDerivatives ReducedReactorDae::derivativesFromVector(
    const numerics::Vector& vector
    ) const
{
    if (vector.size() != DerivativeCount)
    {
        throw std::invalid_argument("Reduced DAE derivative vector has invalid size");
    }

    ReducedReactorDaeDerivatives derivatives;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        derivatives.massDerivativesKgPerS[i] =
            vector[i];
    }

    return derivatives;
}

numerics::Vector ReducedReactorDae::residualToVector(
    const ReducedReactorDaeResidual& residual
    ) const
{
    numerics::Vector vector(ResidualCount, 0.0);

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        vector[i] =
            residual.massResidualsKgPerS[i];
    }

    vector[ComponentCount] =
        residual.volumeResidualM3;

    return vector;
}

numerics::Vector ReducedReactorDae::computeResidualVector(
    const numerics::Vector& variablesVector,
    const numerics::Vector& derivativesVector,
    double temperatureC
    ) const
{
    ReducedReactorDaeVariables variables =
        variablesFromVector(
            variablesVector,
            temperatureC
            );

    ReducedReactorDaeDerivatives derivatives =
        derivativesFromVector(
            derivativesVector
            );

    ReducedReactorDaeResidual residual =
        computeResidual(
            variables,
            derivatives
            );

    return residualToVector(
        residual
        );
}

std::string ReducedReactorDae::residualToString(
    const ReducedReactorDaeResidual& residual
    ) const
{
    std::ostringstream out;

    out << std::fixed << std::setprecision(10);

    const MaterialList& materials =
        residual.massBalance.evaluation.flash.status == FlashStatus::NotConverged
            ? MaterialList{}
            : MaterialList{};

    (void)materials;

    out << "Reduced reactor DAE residual:\n";

    out << "  Mass residuals:\n";

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "    R_M[" << i << "] = "
            << residual.massResidualsKgPerS[i]
            << " kg/s\n";
    }

    out << "  Volume residual:\n";
    out << "    R_V = "
        << residual.volumeResidualM3
        << " m3\n";

    out << "  Underlying mass balance:\n";

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "    f_M[" << i << "] = "
            << residual.massBalance.massDerivativesKgPerS[i]
            << " kg/s\n";
    }

    out << "  Flash state:\n";
    out << "    status = "
        << flashStatusToString(residual.massBalance.evaluation.flash.status)
        << "\n";
    out << "    beta   = "
        << residual.massBalance.evaluation.flash.beta
        << "\n";
    out << "    P      = "
        << residual.massBalance.evaluation.pressurePa / 1.0e5
        << " bar\n";
    out << "    V_calc = "
        << residual.massBalance.evaluation.calculatedVolumeM3
        << " m3\n";

    return out.str();
}

ReducedReactorDaeStepResult ReducedReactorDae::implicitEulerStep(
    const ReactorState& oldState,
    double timeStepS
    ) const
{
    if (timeStepS <= 0.0)
    {
        throw std::invalid_argument("Time step must be positive");
    }

    ReducedReactorDaeVariables oldVariables =
        variablesFromState(oldState);

    numerics::Vector oldVariableVector =
        variablesToVector(oldVariables);

    // Начальное приближение для Newton.
    // Берем старое состояние и делаем грубый явный прогноз по массам.
    MassBalanceResult oldMassBalance =
        reactorModel_.computeMassBalance(oldState);

    ReactorState predictedState =
        oldState;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const Component component =
            static_cast<Component>(i);

        const double predictedMass =
            oldState.massKg(component)
            + timeStepS * oldMassBalance.massDerivativesKgPerS[i];

        predictedState.setMassKg(
            component,
            std::max(0.0, predictedMass)
            );
    }

    // Подбираем давление для прогнозного состояния.
    // Это дает Newton хорошую стартовую точку.
    PressureInitializationResult predictedPressure =
        reactorModel_.initializePressureForVolumeNewton(
            predictedState,
            1.0e4,
            1.0e7,
            1.0e-6,
            15,
            true
            );

    ReactorState initialGuessState =
        predictedPressure.converged
            ? predictedPressure.initializedState
            : oldState;

    ReducedReactorDaeVariables initialGuessVariables =
        variablesFromState(initialGuessState);

    numerics::Vector initialGuess =
        variablesToVector(initialGuessVariables);

    const double temperatureC =
        oldVariables.temperatureC;

    auto residualFunction =
        [this, oldVariableVector, timeStepS, temperatureC]
        (const numerics::Vector& newVariablesVector)
    {
        numerics::Vector derivativesVector(
            DerivativeCount,
            0.0
            );

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            derivativesVector[i] =
                (
                    newVariablesVector[i]
                    - oldVariableVector[i]
                    )
                / timeStepS;
        }

        return computeResidualVector(
            newVariablesVector,
            derivativesVector,
            temperatureC
            );
    };

    numerics::NewtonOptions options;

    options.maxIterations = 30;
    options.residualTolerance = 1.0e-8;
    options.stepTolerance = 1.0e-10;
    options.finiteDifferenceStep = 1.0e-6;
    options.minDampingFactor = 1.0e-6;

    // Ограничения:
    // массы >= 0
    // давление >= 0.01 bar
    options.lowerBounds =
        numerics::Vector{
            0.0,
            0.0,
            0.0,
            0.01
        };

    numerics::NewtonResult newtonResult =
        numerics::NewtonSolver::solve(
            residualFunction,
            initialGuess,
            options
            );

    ReducedReactorDaeStepResult result;

    result.converged =
        newtonResult.converged;

    result.iterations =
        newtonResult.iterations;

    result.timeStepS =
        timeStepS;

    result.variables =
        newtonResult.solution;

    result.residual =
        newtonResult.residual;

    result.residualNorm =
        newtonResult.residualNorm;

    result.message =
        newtonResult.message;

    result.method = "implicit_euler";

    ReducedReactorDaeVariables finalVariables =
        variablesFromVector(
            newtonResult.solution,
            temperatureC
            );

    result.state =
        stateFromVariables(
            finalVariables
            );

    return result;
}

std::string ReducedReactorDae::stepResultToString(
    const ReducedReactorDaeStepResult& result
    ) const
{
    std::ostringstream out;

    out << std::fixed << std::setprecision(10);

    out << "Reduced DAE step result:\n";
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

ReducedReactorDaeStepResult ReducedReactorDae::radauIIA2Step(
    const ReactorState& oldState,
    double timeStepS
    ) const
{
    if (timeStepS <= 0.0)
    {
        throw std::invalid_argument("Time step must be positive");
    }

    constexpr int StageCount = 2;
    constexpr std::size_t StageUnknownCount =
        ComponentCount + 1;

    constexpr double a11 = 5.0 / 12.0;
    constexpr double a12 = -1.0 / 12.0;
    constexpr double a21 = 3.0 / 4.0;
    constexpr double a22 = 1.0 / 4.0;

    constexpr double c1 = 1.0 / 3.0;
    constexpr double c2 = 1.0;

    ReducedReactorDaeVariables oldVariables =
        variablesFromState(oldState);

    numerics::Vector oldVariableVector =
        variablesToVector(oldVariables);

    const double temperatureC =
        oldVariables.temperatureC;

    MassBalanceResult oldMassBalance =
        reactorModel_.computeMassBalance(oldState);

    auto makePressureGuess =
        [this, &oldState](const std::array<double, ComponentCount>& masses)
    {
        ReactorState guessState = oldState;

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            const Component component =
                static_cast<Component>(i);

            guessState.setMassKg(
                component,
                std::max(0.0, masses[i])
                );
        }

        PressureInitializationResult pressureGuess =
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
        const double c =
            stage == 0 ? c1 : c2;

        std::array<double, ComponentCount> stageMassGuess{};

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            stageMassGuess[i] =
                oldVariableVector[i]
                + timeStepS
                      * c
                      * oldMassBalance.massDerivativesKgPerS[i];

            initialGuess[stage * StageUnknownCount + i] =
                oldMassBalance.massDerivativesKgPerS[i];
        }

        initialGuess[stage * StageUnknownCount + ComponentCount] =
            makePressureGuess(stageMassGuess);
    }

    auto stageDerivative =
        [](const numerics::Vector& unknowns, int stage, std::size_t component)
    {
        constexpr std::size_t blockSize =
            ComponentCount + 1;

        return unknowns[
            static_cast<std::size_t>(stage) * blockSize
            + component
        ];
    };

    auto stagePressure =
        [](const numerics::Vector& unknowns, int stage)
    {
        constexpr std::size_t blockSize =
            ComponentCount + 1;

        return unknowns[
            static_cast<std::size_t>(stage) * blockSize
            + ComponentCount
        ];
    };

    auto residualFunction =
        [this,
         oldVariableVector,
         timeStepS,
         temperatureC,
         stageDerivative,
         stagePressure]
        (const numerics::Vector& unknowns)
    {
        constexpr std::size_t totalUnknownCount =
            2 * (ComponentCount + 1);

        if (unknowns.size() != totalUnknownCount)
        {
            throw std::runtime_error("Radau IIA 2 unknown vector has invalid size");
        }

        numerics::Vector totalResidual(
            totalUnknownCount,
            0.0
            );

        for (int stage = 0; stage < 2; ++stage)
        {
            numerics::Vector variablesVector(
                VariableCount,
                0.0
                );

            numerics::Vector derivativesVector(
                DerivativeCount,
                0.0
                );

            for (std::size_t i = 0; i < ComponentCount; ++i)
            {
                const double k1 =
                    stageDerivative(unknowns, 0, i);

                const double k2 =
                    stageDerivative(unknowns, 1, i);

                double mass = oldVariableVector[i];

                if (stage == 0)
                {
                    mass +=
                        timeStepS
                        * (
                            (5.0 / 12.0) * k1
                            + (-1.0 / 12.0) * k2
                            );
                }
                else
                {
                    mass +=
                        timeStepS
                        * (
                            (3.0 / 4.0) * k1
                            + (1.0 / 4.0) * k2
                            );
                }

                variablesVector[i] =
                    mass;

                derivativesVector[i] =
                    stageDerivative(
                        unknowns,
                        stage,
                        i
                        );
            }

            variablesVector[ComponentCount] =
                stagePressure(
                    unknowns,
                    stage
                    );

            numerics::Vector stageResidual =
                computeResidualVector(
                    variablesVector,
                    derivativesVector,
                    temperatureC
                    );

            for (std::size_t i = 0; i < ResidualCount; ++i)
            {
                totalResidual[
                    static_cast<std::size_t>(stage) * ResidualCount
                    + i
                ] = stageResidual[i];
            }
        }

        return totalResidual;
    };

    numerics::NewtonOptions options;

    options.maxIterations = 40;
    options.residualTolerance = 1.0e-8;
    options.stepTolerance = 1.0e-10;
    options.finiteDifferenceStep = 1.0e-6;
    options.minDampingFactor = 1.0e-6;

    const double inf =
        std::numeric_limits<double>::infinity();

    options.lowerBounds =
        numerics::Vector{
            -inf, -inf, -inf, 0.01,
            -inf, -inf, -inf, 0.01
        };

    numerics::NewtonResult newtonResult =
        numerics::NewtonSolver::solve(
            residualFunction,
            initialGuess,
            options
            );

    ReducedReactorDaeStepResult result;

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
        "radau_iia_2";

    // Для Radau IIA последняя стадия совпадает с состоянием в конце шага.
    // Поэтому M_new = old M + h * (3/4*K1 + 1/4*K2)
    // и P_new = P второй стадии.
    ReducedReactorDaeVariables finalVariables;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const double k1 =
            stageDerivative(
                newtonResult.solution,
                0,
                i
                );

        const double k2 =
            stageDerivative(
                newtonResult.solution,
                1,
                i
                );

        finalVariables.massesKg[i] =
            oldVariableVector[i]
            + timeStepS
                  * (
                      (3.0 / 4.0) * k1
                      + (1.0 / 4.0) * k2
                      );
    }

    finalVariables.pressureBar =
        stagePressure(
            newtonResult.solution,
            1
            );

    finalVariables.temperatureC =
        temperatureC;

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

ReducedReactorDaeStepResult ReducedReactorDae::radauIIA3Step(
    const ReactorState& oldState,
    double timeStepS
    ) const
{
    if (timeStepS <= 0.0)
    {
        throw std::invalid_argument("Time step must be positive");
    }

    constexpr int StageCount = 3;
    constexpr std::size_t StageUnknownCount =
        ComponentCount + 1;

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

    ReducedReactorDaeVariables oldVariables =
        variablesFromState(oldState);

    numerics::Vector oldVariableVector =
        variablesToVector(oldVariables);

    const double temperatureC =
        oldVariables.temperatureC;

    MassBalanceResult oldMassBalance =
        reactorModel_.computeMassBalance(oldState);

    auto makePressureGuess =
        [this, &oldState](const std::array<double, ComponentCount>& masses)
    {
        ReactorState guessState = oldState;

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            const Component component =
                static_cast<Component>(i);

            guessState.setMassKg(
                component,
                std::max(0.0, masses[i])
                );
        }

        PressureInitializationResult pressureGuess =
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
                oldVariableVector[i]
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
            makePressureGuess(stageMassGuess);
    }

    auto stageDerivative =
        [](const numerics::Vector& unknowns, int stage, std::size_t component)
    {
        constexpr std::size_t blockSize =
            ComponentCount + 1;

        return unknowns[
            static_cast<std::size_t>(stage) * blockSize
            + component
        ];
    };

    auto stagePressure =
        [](const numerics::Vector& unknowns, int stage)
    {
        constexpr std::size_t blockSize =
            ComponentCount + 1;

        return unknowns[
            static_cast<std::size_t>(stage) * blockSize
            + ComponentCount
        ];
    };

    auto residualFunction =
        [this,
         oldVariableVector,
         timeStepS,
         temperatureC,
         stageDerivative,
         stagePressure,
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
            3 * (ComponentCount + 1);

        if (unknowns.size() != totalUnknownCount)
        {
            throw std::runtime_error("Radau IIA 3 unknown vector has invalid size");
        }

        numerics::Vector totalResidual(
            totalUnknownCount,
            0.0
            );

        for (int stage = 0; stage < 3; ++stage)
        {
            numerics::Vector variablesVector(
                VariableCount,
                0.0
                );

            numerics::Vector derivativesVector(
                DerivativeCount,
                0.0
                );

            for (std::size_t i = 0; i < ComponentCount; ++i)
            {
                const double k1 =
                    stageDerivative(unknowns, 0, i);

                const double k2 =
                    stageDerivative(unknowns, 1, i);

                const double k3 =
                    stageDerivative(unknowns, 2, i);

                double mass =
                    oldVariableVector[i];

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

                variablesVector[i] =
                    mass;

                derivativesVector[i] =
                    stageDerivative(
                        unknowns,
                        stage,
                        i
                        );
            }

            variablesVector[ComponentCount] =
                stagePressure(
                    unknowns,
                    stage
                    );

            numerics::Vector stageResidual =
                computeResidualVector(
                    variablesVector,
                    derivativesVector,
                    temperatureC
                    );

            for (std::size_t i = 0; i < ResidualCount; ++i)
            {
                totalResidual[
                    static_cast<std::size_t>(stage) * ResidualCount
                    + i
                ] =
                    stageResidual[i];
            }
        }

        return totalResidual;
    };

    numerics::NewtonOptions options;

    options.maxIterations = 50;
    options.residualTolerance = 1.0e-8;
    options.stepTolerance = 1.0e-10;
    options.finiteDifferenceStep = 1.0e-6;
    options.minDampingFactor = 1.0e-6;

    const double inf =
        std::numeric_limits<double>::infinity();

    options.lowerBounds =
        numerics::Vector{
            -inf, -inf, -inf, 0.01,
            -inf, -inf, -inf, 0.01,
            -inf, -inf, -inf, 0.01
        };

    numerics::NewtonResult newtonResult =
        numerics::NewtonSolver::solve(
            residualFunction,
            initialGuess,
            options
            );

    ReducedReactorDaeStepResult result;

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
        "radau_iia_3";

    ReducedReactorDaeVariables finalVariables;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const double k1 =
            stageDerivative(
                newtonResult.solution,
                0,
                i
                );

        const double k2 =
            stageDerivative(
                newtonResult.solution,
                1,
                i
                );

        const double k3 =
            stageDerivative(
                newtonResult.solution,
                2,
                i
                );

        finalVariables.massesKg[i] =
            oldVariableVector[i]
            + timeStepS
                  * (
                      a31 * k1
                      + a32 * k2
                      + a33 * k3
                      );
    }

    finalVariables.pressureBar =
        stagePressure(
            newtonResult.solution,
            2
            );

    finalVariables.temperatureC =
        temperatureC;

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
