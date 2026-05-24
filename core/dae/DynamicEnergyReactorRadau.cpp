#include "DynamicEnergyReactorRadau.h"

#include <array>
#include <cmath>
#include <stdexcept>

namespace {

constexpr int StageCount = 3;

constexpr double Sqrt6 = 2.4494897427831780982;

// Radau IIA 3-stage, order 5
const std::array<double, StageCount> C = {
    (4.0 - Sqrt6) / 10.0,
    (4.0 + Sqrt6) / 10.0,
    1.0
};

const std::array<std::array<double, StageCount>, StageCount> A = {{
    {
        (88.0 - 7.0 * Sqrt6) / 360.0,
        (296.0 - 169.0 * Sqrt6) / 1800.0,
        (-2.0 + 3.0 * Sqrt6) / 225.0
    },
    {
        (296.0 + 169.0 * Sqrt6) / 1800.0,
        (88.0 + 7.0 * Sqrt6) / 360.0,
        (-2.0 - 3.0 * Sqrt6) / 225.0
    },
    {
        (16.0 - Sqrt6) / 36.0,
        (16.0 + Sqrt6) / 36.0,
        1.0 / 9.0
    }
}};

const std::array<double, StageCount> B = {
    (16.0 - Sqrt6) / 36.0,
    (16.0 + Sqrt6) / 36.0,
    1.0 / 9.0
};

} // namespace

DynamicEnergyReactorRadau::DynamicEnergyReactorRadau(
    const DynamicEnergyReactorDae& dae,
    const DynamicWellStirredReactorModel& reactorModel
    )
    : dae_(dae),
    reactorModel_(reactorModel)
{
}

std::size_t DynamicEnergyReactorRadau::differentialCount() const {
    return dae_.layout().componentCount() + 1;
}

std::size_t DynamicEnergyReactorRadau::differentialVariableIndex(
    std::size_t differentialIndex
    ) const {
    const DynamicStateLayout& layout =
        dae_.layout();

    const std::size_t n =
        layout.componentCount();

    if (differentialIndex < n) {
        return layout.mass(differentialIndex);
    }

    if (differentialIndex == n) {
        return layout.energy();
    }

    throw std::runtime_error(
        "DynamicEnergyReactorRadau: differential index out of range"
        );
}

std::size_t DynamicEnergyReactorRadau::stageBlockSize() const {
    return dae_.layout().variableCount() + differentialCount();
}

std::size_t DynamicEnergyReactorRadau::stageOffset(
    std::size_t stageIndex
    ) const {
    return stageIndex * stageBlockSize();
}

std::size_t DynamicEnergyReactorRadau::stageVariableIndex(
    std::size_t stageIndex,
    std::size_t variableIndex
    ) const {
    return stageOffset(stageIndex) + variableIndex;
}

std::size_t DynamicEnergyReactorRadau::stageDerivativeIndex(
    std::size_t stageIndex,
    std::size_t differentialIndex
    ) const {
    return stageOffset(stageIndex)
    + dae_.layout().variableCount()
        + differentialIndex;
}

DynamicVector DynamicEnergyReactorRadau::extractStageVariables(
    const DynamicVector& unknowns,
    std::size_t stageIndex
    ) const {
    const std::size_t variableCount =
        dae_.layout().variableCount();

    DynamicVector variables(variableCount, 0.0);

    for (std::size_t i = 0; i < variableCount; ++i) {
        variables[i] =
            unknowns[stageVariableIndex(stageIndex, i)];
    }

    return variables;
}

DynamicVector DynamicEnergyReactorRadau::extractStageDerivativesFull(
    const DynamicVector& unknowns,
    std::size_t stageIndex
    ) const {
    const DynamicStateLayout& layout = dae_.layout();

    DynamicVector derivatives(
        layout.variableCount(),
        0.0
        );

    for (std::size_t j = 0; j < differentialCount(); ++j) {
        const std::size_t variableIndex =
            differentialVariableIndex(j);

        derivatives[variableIndex] =
            unknowns[stageDerivativeIndex(stageIndex, j)];
    }

    return derivatives;
}

DynamicVector DynamicEnergyReactorRadau::buildInitialGuess(
    const DynamicVector& oldVariables,
    double timeStepS
    ) const {
    const DynamicStateLayout& layout = dae_.layout();

    const std::size_t variableCount =
        layout.variableCount();

    const std::size_t blockSize =
        stageBlockSize();

    DynamicVector guess(
        StageCount * blockSize,
        0.0
        );

    const DynamicEnergyReactorDaeEvaluation oldEvaluation =
        dae_.evaluate(oldVariables);

    const DynamicMassBalanceResult oldMassBalance =
        reactorModel_.computeMassBalance(
            oldEvaluation.pressurePa,
            oldEvaluation.massesKg,
            oldEvaluation.gibbsOde.yVapor
            );

    const DynamicEnergyBalanceResult oldEnergyBalance =
        reactorModel_.computeEnergyBalance(
            oldEvaluation.pressurePa,
            oldEvaluation.temperatureK,
            oldMassBalance.outletMolarFlowKmolPerS,
            oldEvaluation.gibbsOde.yVapor
            );

    DynamicVector oldDerivativeFull(
        variableCount,
        0.0
        );

    for (std::size_t i = 0; i < layout.componentCount(); ++i) {
        oldDerivativeFull[layout.mass(i)] =
            oldMassBalance.massRatesKgPerS[i];
    }

    oldDerivativeFull[layout.energy()] =
        oldEnergyBalance.energyRateW;


    for (std::size_t stage = 0; stage < StageCount; ++stage) {
        for (std::size_t i = 0; i < variableCount; ++i) {
            guess[stageVariableIndex(stage, i)] =
                oldVariables[i];
        }

        // Начальное приближение для дифференциальных переменных:
        // y_stage ≈ y_old + c_i * dt * f_old
        for (std::size_t j = 0; j < differentialCount(); ++j) {
            const std::size_t variableIndex =
                differentialVariableIndex(j);

            guess[stageVariableIndex(stage, variableIndex)] =
                oldVariables[variableIndex]
                + C[stage] * timeStepS
                      * oldDerivativeFull[variableIndex];

            guess[stageDerivativeIndex(stage, j)] =
                oldDerivativeFull[variableIndex];
        }
    }

    return guess;
}

DynamicVector DynamicEnergyReactorRadau::buildResidual(
    const DynamicVector& unknowns,
    const DynamicVector& oldVariables,
    double timeStepS
    ) const {
    const DynamicStateLayout& layout = dae_.layout();

    const std::size_t variableCount =
        layout.variableCount();

    const std::size_t dCount =
        differentialCount();

    const std::size_t blockSize =
        stageBlockSize();

    DynamicVector residual(
        StageCount * blockSize,
        0.0
        );

    for (std::size_t stage = 0; stage < StageCount; ++stage) {
        const DynamicVector stageVariables =
            extractStageVariables(
                unknowns,
                stage
                );

        const DynamicVector stageDerivatives =
            extractStageDerivativesFull(
                unknowns,
                stage
                );

        const DynamicVector daeResidual =
            dae_.computeCoupledResidualVector(
                stageVariables,
                stageDerivatives,
                reactorModel_
                );

        for (std::size_t i = 0; i < variableCount; ++i) {
            residual[stageOffset(stage) + i] =
                daeResidual[i];
        }

        for (std::size_t j = 0; j < dCount; ++j) {
            const std::size_t variableIndex =
                differentialVariableIndex(j);

            double collocationValue =
                oldVariables[variableIndex];

            for (std::size_t column = 0; column < StageCount; ++column) {
                collocationValue +=
                    timeStepS
                    * A[stage][column]
                    * unknowns[stageDerivativeIndex(column, j)];
            }

            residual[
                stageOffset(stage)
                + variableCount
                + j
            ] =
                stageVariables[variableIndex]
                - collocationValue;
        }
    }

    return residual;
}

DynamicVector DynamicEnergyReactorRadau::buildNewVariables(
    const DynamicVector& oldVariables,
    const DynamicVector& solvedUnknowns,
    double timeStepS
    ) const {
    const DynamicStateLayout& layout = dae_.layout();

    DynamicVector newVariables =
        oldVariables;

    for (std::size_t j = 0; j < differentialCount(); ++j) {
        const std::size_t variableIndex =
            differentialVariableIndex(j);

        double value =
            oldVariables[variableIndex];

        for (std::size_t stage = 0; stage < StageCount; ++stage) {
            value +=
                timeStepS
                * B[stage]
                * solvedUnknowns[
                    stageDerivativeIndex(stage, j)
            ];
        }

        newVariables[variableIndex] =
            value;
    }

    // У Radau IIA последняя стадия находится в конце шага: c_3 = 1.
    // Поэтому алгебраические переменные берем из последней стадии.
    const std::size_t lastStage =
        StageCount - 1;

    newVariables[layout.pressure()] =
        solvedUnknowns[
            stageVariableIndex(lastStage, layout.pressure())
    ];

    newVariables[layout.temperature()] =
        solvedUnknowns[
            stageVariableIndex(lastStage, layout.temperature())
    ];

    for (std::size_t i = 0; i < layout.componentCount(); ++i) {
        newVariables[layout.liquidSplit(i)] =
            solvedUnknowns[
                stageVariableIndex(
                    lastStage,
                    layout.liquidSplit(i)
                    )
        ];
    }

    newVariables[layout.liquidMolarVolume()] =
        solvedUnknowns[
            stageVariableIndex(
                lastStage,
                layout.liquidMolarVolume()
                )
    ];

    newVariables[layout.vaporMolarVolume()] =
        solvedUnknowns[
            stageVariableIndex(
                lastStage,
                layout.vaporMolarVolume()
                )
    ];



    constexpr double NegativeMassToleranceKg = 1.0e-8;

    for (std::size_t i = 0; i < layout.componentCount(); ++i) {
        const std::size_t massIndex =
            layout.mass(i);

        if (newVariables[massIndex] < -NegativeMassToleranceKg) {
            throw std::runtime_error(
                "DynamicEnergyReactorRadau: negative component mass after step"
                );
        }

        if (newVariables[massIndex] < 0.0) {
            newVariables[massIndex] = 0.0;
        }
    }

    for (std::size_t i = 0; i < layout.componentCount(); ++i) {
        const std::size_t splitIndex =
            layout.liquidSplit(i);

        newVariables[splitIndex] =
            std::max(
                1.0e-12,
                std::min(
                    newVariables[splitIndex],
                    1.0 - 1.0e-12
                    )
                );
    }

    newVariables[layout.liquidMolarVolume()] =
        std::max(
            newVariables[layout.liquidMolarVolume()],
            1.0e-9
            );

    newVariables[layout.vaporMolarVolume()] =
        std::max(
            newVariables[layout.vaporMolarVolume()],
            1.0e-9
            );

    return newVariables;
}

DynamicRadauStepResult DynamicEnergyReactorRadau::radauIIA3Step(
    const DynamicVector& oldVariables,
    double timeStepS,
    const DynamicRadauOptions& options
    ) const {
    if (timeStepS <= 0.0) {
        throw std::runtime_error(
            "DynamicEnergyReactorRadau: time step must be positive"
            );
    }

    const std::size_t variableCount =
        dae_.layout().variableCount();

    if (oldVariables.size() != variableCount) {
        throw std::runtime_error(
            "DynamicEnergyReactorRadau: invalid oldVariables size"
            );
    }

    DynamicRadauStepResult stepResult;
    stepResult.oldVariables = oldVariables;

    numerics::NewtonOptions newtonOptions;
    newtonOptions.residualTolerance =
        options.residualTolerance;
    newtonOptions.stepTolerance =
        options.stepTolerance;
    newtonOptions.maxIterations =
        options.maxNewtonIterations;

    const DynamicVector initialGuess =
        buildInitialGuess(
            oldVariables,
            timeStepS
            );

    auto residualFunction =
        [&](const numerics::Vector& unknowns) {
            return buildResidual(
                unknowns,
                oldVariables,
                timeStepS
                );
        };

    const numerics::NewtonResult newtonResult =
        numerics::NewtonSolver::solve(
            residualFunction,
            initialGuess,
            newtonOptions
            );

    stepResult.converged =
        newtonResult.converged;

    stepResult.iterations =
        newtonResult.iterations;

    stepResult.residualNorm =
        newtonResult.residualNorm;

    stepResult.message =
        newtonResult.message;

    if (!newtonResult.converged) {
        stepResult.newVariables =
            oldVariables;
        return stepResult;
    }

    stepResult.newVariables =
        buildNewVariables(
            oldVariables,
            newtonResult.solution,
            timeStepS
            );

    return stepResult;
}
