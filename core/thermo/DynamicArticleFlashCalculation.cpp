#include "DynamicArticleFlashCalculation.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace {

constexpr double GasConstantJPerKmolK = 8314.46261815324;

}

DynamicArticleFlashCalculation::DynamicArticleFlashCalculation(
    const DynamicPengRobinsonEOS& eos
    )
    : eos_(eos),
    odeSystem_(eos)
{
}

double DynamicArticleFlashCalculation::clampValue(
    double value,
    double minValue,
    double maxValue
    ) const {
    return std::max(
        minValue,
        std::min(value, maxValue)
        );
}

DynamicComposition DynamicArticleFlashCalculation::normalizeComposition(
    const DynamicComposition& composition
    ) const {
    if (composition.size() != eos_.componentCount()) {
        throw std::runtime_error(
            "DynamicArticleFlashCalculation: composition size mismatch"
            );
    }

    double sum = 0.0;

    for (double value : composition) {
        if (value < 0.0) {
            throw std::runtime_error(
                "DynamicArticleFlashCalculation: negative composition value"
                );
        }

        sum += value;
    }

    if (sum <= 0.0) {
        throw std::runtime_error(
            "DynamicArticleFlashCalculation: composition sum must be positive"
            );
    }

    DynamicComposition z =
        composition;

    for (double& value : z) {
        value /= sum;
    }

    return z;
}

DynamicComposition DynamicArticleFlashCalculation::initialLiquidSplit(
    const DynamicComposition& z,
    const DynamicArticleFlashOptions& options
    ) const {
    DynamicComposition split(
        z.size(),
        clampValue(
            options.initialLiquidSplit,
            options.minSplit,
            1.0 - options.minSplit
            )
        );

    for (std::size_t i = 0; i < z.size(); ++i) {
        if (z[i] <= options.minSplit) {
            split[i] = options.minSplit;
        }
    }

    return split;
}

void DynamicArticleFlashCalculation::initialMolarVolumes(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& z,
    double& liquidMolarVolumeM3PerKmol,
    double& vaporMolarVolumeM3PerKmol
    ) const {
    const DynamicZFactors roots =
        eos_.computeZFactors(
            pressurePa,
            temperatureK,
            z
            );

    double liquidZ = roots.liquid;
    double vaporZ = roots.vapor;

    if (roots.roots.size() == 1) {
        const double singleZ =
            roots.roots.front();

        if (singleZ < 0.3) {
            liquidZ = singleZ;
            vaporZ = std::max(1.0, singleZ * 10.0);
        } else {
            liquidZ = std::max(0.02, singleZ * 0.1);
            vaporZ = singleZ;
        }
    }

    liquidMolarVolumeM3PerKmol =
        liquidZ
        * GasConstantJPerKmolK
        * temperatureK
        / pressurePa;

    vaporMolarVolumeM3PerKmol =
        vaporZ
        * GasConstantJPerKmolK
        * temperatureK
        / pressurePa;

    liquidMolarVolumeM3PerKmol =
        std::max(
            liquidMolarVolumeM3PerKmol,
            1.0e-9
            );

    vaporMolarVolumeM3PerKmol =
        std::max(
            vaporMolarVolumeM3PerKmol,
            1.0e-9
            );
}

bool DynamicArticleFlashCalculation::isSingleLiquidBoundary(
    const DynamicGibbsOdeState& state,
    const DynamicArticleFlashOptions& options
    ) const {
    return state.betaVapor <= options.singlePhaseTolerance;
}

bool DynamicArticleFlashCalculation::isSingleVaporBoundary(
    const DynamicGibbsOdeState& state,
    const DynamicArticleFlashOptions& options
    ) const {
    return state.betaLiquid <= options.singlePhaseTolerance;
}

bool DynamicArticleFlashCalculation::isConvergedTwoPhase(
    const DynamicGibbsOdeState& state,
    double pressurePa,
    const DynamicArticleFlashOptions& options,
    double& maxPressureResidual,
    double& maxRateResidual
    ) const {
    const double liquidPressureResidual =
        std::abs(
            state.liquidPressureFromEosPa - pressurePa
            ) / pressurePa;

    const double vaporPressureResidual =
        std::abs(
            state.vaporPressureFromEosPa - pressurePa
            ) / pressurePa;

    maxPressureResidual =
        std::max(
            liquidPressureResidual,
            vaporPressureResidual
            );

    maxRateResidual =
        std::max(
            std::abs(state.liquidMolarVolumeRate),
            std::abs(state.vaporMolarVolumeRate)
            );

    for (double rate : state.liquidSplitRates) {
        maxRateResidual =
            std::max(
                maxRateResidual,
                std::abs(rate)
                );
    }

    return
        state.maxFugacityResidual <= options.fugacityTolerance
        &&
        maxPressureResidual <= options.pressureTolerance;
}

DynamicArticleFlashResult DynamicArticleFlashCalculation::buildResult(
    double pressurePa,
    double temperatureK,
    const DynamicGibbsOdeState& state,
    int iterations,
    bool converged,
    const std::string& status,
    const DynamicArticleFlashOptions& options
    ) const {
    DynamicArticleFlashResult result;

    result.converged = converged;
    result.status = status;
    result.iterations = iterations;

    result.betaLiquid = state.betaLiquid;
    result.betaVapor = state.betaVapor;

    result.liquidSplit = state.liquidSplit;

    result.xLiquid = state.xLiquid;
    result.yVapor = state.yVapor;

    result.phiLiquid = state.phiLiquid;
    result.phiVapor = state.phiVapor;

    result.fugacityLogResiduals =
        state.fugacityLogResiduals;

    result.liquidMolarVolumeM3PerKmol =
        state.liquidMolarVolumeM3PerKmol;

    result.vaporMolarVolumeM3PerKmol =
        state.vaporMolarVolumeM3PerKmol;

    result.liquidZ = state.liquidZ;
    result.vaporZ = state.vaporZ;

    result.liquidPressureFromEosPa =
        state.liquidPressureFromEosPa;

    result.vaporPressureFromEosPa =
        state.vaporPressureFromEosPa;

    result.gibbsDeviationDimensionless =
        state.gibbsDeviationDimensionless;

    result.maxFugacityResidual =
        state.maxFugacityResidual;

    const double liquidPressureResidual =
        std::abs(
            state.liquidPressureFromEosPa - pressurePa
            ) / pressurePa;

    const double vaporPressureResidual =
        std::abs(
            state.vaporPressureFromEosPa - pressurePa
            ) / pressurePa;

    result.maxPressureResidual =
        std::max(
            liquidPressureResidual,
            vaporPressureResidual
            );

    result.maxRateResidual =
        std::max(
            std::abs(state.liquidMolarVolumeRate),
            std::abs(state.vaporMolarVolumeRate)
            );

    for (double rate : state.liquidSplitRates) {
        result.maxRateResidual =
            std::max(
                result.maxRateResidual,
                std::abs(rate)
                );
    }

    if (isSingleLiquidBoundary(state, options)) {
        result.phaseState =
            DynamicArticlePhaseState::SingleLiquid;

        result.betaLiquid = 1.0;
        result.betaVapor = 0.0;
        result.status =
            converged
                ? "Single liquid by Gibbs ODE boundary minimum"
                : status;
    } else if (isSingleVaporBoundary(state, options)) {
        result.phaseState =
            DynamicArticlePhaseState::SingleVapor;

        result.betaLiquid = 0.0;
        result.betaVapor = 1.0;
        result.status =
            converged
                ? "Single vapor by Gibbs ODE boundary minimum"
                : status;
    } else if (converged) {
        result.phaseState =
            DynamicArticlePhaseState::TwoPhase;
    } else {
        result.phaseState =
            DynamicArticlePhaseState::NotConverged;
    }

    return result;
}

DynamicArticleFlashResult DynamicArticleFlashCalculation::calculate(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& overallComposition,
    const DynamicArticleFlashOptions& options
    ) const {
    if (pressurePa <= 0.0 || temperatureK <= 0.0) {
        throw std::runtime_error(
            "DynamicArticleFlashCalculation: pressure and temperature must be positive"
            );
    }

    if (options.pseudoTimeStep <= 0.0) {
        throw std::runtime_error(
            "DynamicArticleFlashCalculation: pseudo time step must be positive"
            );
    }

    const DynamicComposition z =
        normalizeComposition(overallComposition);

    DynamicComposition split =
        initialLiquidSplit(
            z,
            options
            );

    double liquidVolume = 0.0;
    double vaporVolume = 0.0;

    initialMolarVolumes(
        pressurePa,
        temperatureK,
        z,
        liquidVolume,
        vaporVolume
        );

    DynamicGibbsOdeOptions odeOptions;
    odeOptions.minSplit = options.minSplit;
    odeOptions.compositionMobility =
        options.compositionMobility;
    odeOptions.volumeMobility =
        options.volumeMobility;

    DynamicGibbsOdeState state;

    for (int iteration = 0;
         iteration < options.maxIterations;
         ++iteration) {

        state =
            odeSystem_.evaluate(
                pressurePa,
                temperatureK,
                1.0,
                z,
                split,
                liquidVolume,
                vaporVolume,
                odeOptions
                );

        double maxPressureResidual = 0.0;
        double maxRateResidual = 0.0;

        if (isSingleLiquidBoundary(state, options)) {
            return buildResult(
                pressurePa,
                temperatureK,
                state,
                iteration,
                true,
                "Single liquid by Gibbs ODE boundary minimum",
                options
                );
        }

        if (isSingleVaporBoundary(state, options)) {
            return buildResult(
                pressurePa,
                temperatureK,
                state,
                iteration,
                true,
                "Single vapor by Gibbs ODE boundary minimum",
                options
                );
        }

        if (isConvergedTwoPhase(
                state,
                pressurePa,
                options,
                maxPressureResidual,
                maxRateResidual
                )) {
            return buildResult(
                pressurePa,
                temperatureK,
                state,
                iteration,
                true,
                "Two-phase Gibbs ODE flash converged",
                options
                );
        }

        for (std::size_t i = 0; i < split.size(); ++i) {
            if (z[i] <= options.minSplit) {
                split[i] = options.minSplit;
                continue;
            }

            split[i] +=
                options.pseudoTimeStep
                * state.liquidSplitRates[i];

            split[i] =
                clampValue(
                    split[i],
                    options.minSplit,
                    1.0 - options.minSplit
                    );
        }

        liquidVolume +=
            options.pseudoTimeStep
            * state.liquidMolarVolumeRate;

        vaporVolume +=
            options.pseudoTimeStep
            * state.vaporMolarVolumeRate;

        liquidVolume =
            std::max(
                liquidVolume,
                1.0e-9
                );

        vaporVolume =
            std::max(
                vaporVolume,
                1.0e-9
                );
    }

    state =
        odeSystem_.evaluate(
            pressurePa,
            temperatureK,
            1.0,
            z,
            split,
            liquidVolume,
            vaporVolume,
            odeOptions
            );

    return buildResult(
        pressurePa,
        temperatureK,
        state,
        options.maxIterations,
        false,
        "Gibbs ODE flash did not converge",
        options
        );
}
