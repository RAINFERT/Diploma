#include "DynamicGibbsOdeSystem.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace {

constexpr double GasConstantJPerKmolK = 8314.46261815324;

} // namespace

DynamicGibbsOdeSystem::DynamicGibbsOdeSystem(
    const DynamicPengRobinsonEOS& eos
    )
    : eos_(eos)
{
}

double DynamicGibbsOdeSystem::safePositive(
    double value
    ) const {
    return std::max(value, 1.0e-300);
}

double DynamicGibbsOdeSystem::clampValue(
    double value,
    double minValue,
    double maxValue
    ) const {
    return std::max(minValue, std::min(value, maxValue));
}

DynamicComposition DynamicGibbsOdeSystem::normalizeComposition(
    const DynamicComposition& composition
    ) const {
    if (composition.size() != eos_.componentCount()) {
        throw std::runtime_error(
            "DynamicGibbsOdeSystem: composition size mismatch"
            );
    }

    double sum = 0.0;

    for (double value : composition) {
        if (value < 0.0) {
            throw std::runtime_error(
                "DynamicGibbsOdeSystem: negative composition value"
                );
        }

        sum += value;
    }

    if (sum <= 0.0) {
        throw std::runtime_error(
            "DynamicGibbsOdeSystem: composition sum must be positive"
            );
    }

    DynamicComposition result =
        composition;

    for (double& value : result) {
        value /= sum;
    }

    return result;
}

double DynamicGibbsOdeSystem::pseudoCriticalVolumeM3PerKmol(
    const DynamicComposition& composition
    ) const {
    const DynamicComposition z =
        normalizeComposition(composition);

    double value = 0.0;

    for (std::size_t i = 0; i < z.size(); ++i) {
        const ComponentProperties& component =
            eos_.component(i);

        double vc =
            component.criticalVolumeM3PerKmol;

        if (vc <= 0.0) {
            vc =
                component.criticalCompressibility
                * GasConstantJPerKmolK
                * component.criticalTemperatureK
                / component.criticalPressurePa;
        }

        if (vc <= 0.0) {
            throw std::runtime_error(
                "DynamicGibbsOdeSystem: invalid critical volume"
                );
        }

        value +=
            z[i] * vc;
    }

    return value;
}

DynamicGibbsOdeState DynamicGibbsOdeSystem::evaluate(
    double pressurePa,
    double temperatureK,
    double totalMolesKmol,
    const DynamicComposition& overallComposition,
    const DynamicComposition& liquidSplit,
    double liquidMolarVolumeM3PerKmol,
    double vaporMolarVolumeM3PerKmol,
    const DynamicGibbsOdeOptions& options
    ) const {
    if (pressurePa <= 0.0 || temperatureK <= 0.0) {
        throw std::runtime_error(
            "DynamicGibbsOdeSystem: pressure and temperature must be positive"
            );
    }

    if (totalMolesKmol <= 0.0) {
        throw std::runtime_error(
            "DynamicGibbsOdeSystem: total moles must be positive"
            );
    }

    if (liquidSplit.size() != eos_.componentCount()) {
        throw std::runtime_error(
            "DynamicGibbsOdeSystem: liquidSplit size mismatch"
            );
    }

    if (liquidMolarVolumeM3PerKmol <= 0.0 ||
        vaporMolarVolumeM3PerKmol <= 0.0) {
        throw std::runtime_error(
            "DynamicGibbsOdeSystem: molar volumes must be positive"
            );
    }

    const DynamicComposition z =
        normalizeComposition(overallComposition);

    const std::size_t n =
        eos_.componentCount();

    DynamicGibbsOdeState state;

    state.overallComposition = z;
    state.liquidSplit.assign(n, 0.0);
    state.liquidMolesFraction.assign(n, 0.0);
    state.vaporMolesFraction.assign(n, 0.0);
    state.xLiquid.assign(n, 0.0);
    state.yVapor.assign(n, 0.0);
    state.fugacityLogResiduals.assign(n, 0.0);
    state.liquidSplitRates.assign(n, 0.0);

    state.liquidMolarVolumeM3PerKmol =
        liquidMolarVolumeM3PerKmol;

    state.vaporMolarVolumeM3PerKmol =
        vaporMolarVolumeM3PerKmol;

    for (std::size_t i = 0; i < n; ++i) {
        const double zi =
            z[i];

        if (zi <= 0.0) {
            state.liquidSplit[i] =
                options.minSplit;

            state.liquidMolesFraction[i] = 0.0;
            state.vaporMolesFraction[i] = 0.0;
            continue;
        }

        const double si =
            clampValue(
                liquidSplit[i],
                options.minSplit,
                1.0 - options.minSplit
                );

        state.liquidSplit[i] =
            si;

        state.liquidMolesFraction[i] =
            zi * si;

        state.vaporMolesFraction[i] =
            zi * (1.0 - si);

        state.betaLiquid +=
            state.liquidMolesFraction[i];

        state.betaVapor +=
            state.vaporMolesFraction[i];
    }

    const double betaLiquidSafe =
        safePositive(state.betaLiquid);

    const double betaVaporSafe =
        safePositive(state.betaVapor);

    for (std::size_t i = 0; i < n; ++i) {
        state.xLiquid[i] =
            state.liquidMolesFraction[i] / betaLiquidSafe;

        state.yVapor[i] =
            state.vaporMolesFraction[i] / betaVaporSafe;
    }

    state.liquidZ =
        eos_.computeZFactorFromMolarVolume(
            pressurePa,
            temperatureK,
            liquidMolarVolumeM3PerKmol
            );

    state.vaporZ =
        eos_.computeZFactorFromMolarVolume(
            pressurePa,
            temperatureK,
            vaporMolarVolumeM3PerKmol
            );

    state.phiLiquid =
        eos_.computeFugacityCoefficients(
            pressurePa,
            temperatureK,
            state.xLiquid,
            state.liquidZ
            );

    state.phiVapor =
        eos_.computeFugacityCoefficients(
            pressurePa,
            temperatureK,
            state.yVapor,
            state.vaporZ
            );

    state.liquidPressureFromEosPa =
        eos_.computePressureFromMolarVolume(
            temperatureK,
            state.xLiquid,
            liquidMolarVolumeM3PerKmol
            );

    state.vaporPressureFromEosPa =
        eos_.computePressureFromMolarVolume(
            temperatureK,
            state.yVapor,
            vaporMolarVolumeM3PerKmol
            );

    double maxResidual = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        if (z[i] <= 0.0) {
            state.fugacityLogResiduals[i] = 0.0;
            state.liquidSplitRates[i] = 0.0;
            continue;
        }

        const double fLiquid =
            safePositive(
                state.xLiquid[i]
                * state.phiLiquid[i]
                * pressurePa
                );

        const double fVapor =
            safePositive(
                state.yVapor[i]
                * state.phiVapor[i]
                * pressurePa
                );

        const double residual =
            std::log(fLiquid) - std::log(fVapor);

        state.fugacityLogResiduals[i] =
            residual;

        maxResidual =
            std::max(
                maxResidual,
                std::abs(residual)
                );
    }

    state.maxFugacityResidual =
        maxResidual;

    double g = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        if (z[i] <= 0.0) {
            continue;
        }

        const double nL =
            state.liquidMolesFraction[i];

        const double nV =
            state.vaporMolesFraction[i];

        const double liquidArgument =
            safePositive(
                state.xLiquid[i]
                * state.phiLiquid[i]
                / z[i]
                );

        const double vaporArgument =
            safePositive(
                state.yVapor[i]
                * state.phiVapor[i]
                / z[i]
                );

        g +=
            nL * std::log(liquidArgument)
            +
            nV * std::log(vaporArgument);
    }

    state.gibbsDeviationDimensionless =
        g;

    const double probabilityFactor =
        std::exp(
            -clampValue(
                g,
                -options.maxExpArgument,
                options.maxExpArgument
                )
            );

    for (std::size_t i = 0; i < n; ++i) {
        if (z[i] <= 0.0) {
            state.liquidSplitRates[i] = 0.0;
            continue;
        }

        const double si =
            state.liquidSplit[i];

        const double boundaryFactor =
            si * (1.0 - si);

        state.liquidSplitRates[i] =
            -options.compositionMobility
            * probabilityFactor
            * boundaryFactor
            * state.fugacityLogResiduals[i];
    }

    const double vcLiquid =
        pseudoCriticalVolumeM3PerKmol(
            state.xLiquid
            );

    const double vcVapor =
        pseudoCriticalVolumeM3PerKmol(
            state.yVapor
            );

    const double liquidPressureResidual =
        clampValue(
            (state.liquidPressureFromEosPa - pressurePa) / pressurePa,
            -options.maxPressureResidual,
            options.maxPressureResidual
            );

    const double vaporPressureResidual =
        clampValue(
            (state.vaporPressureFromEosPa - pressurePa) / pressurePa,
            -options.maxPressureResidual,
            options.maxPressureResidual
            );

    state.liquidMolarVolumeRate =
        options.volumeMobility
        * probabilityFactor
        * vcLiquid
        * liquidPressureResidual;

    state.vaporMolarVolumeRate =
        options.volumeMobility
        * probabilityFactor
        * vcVapor
        * vaporPressureResidual;

    return state;
}
