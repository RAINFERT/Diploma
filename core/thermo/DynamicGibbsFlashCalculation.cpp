#include "DynamicGibbsFlashCalculation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace {

double clampValue(
    double value,
    double minValue,
    double maxValue
    ) {
    return std::max(minValue, std::min(value, maxValue));
}

double safePositive(
    double value
    ) {
    constexpr double MinValue = 1.0e-300;
    return std::max(value, MinValue);
}

double maxAbsValue(
    const DynamicComposition& values
    ) {
    double result = 0.0;

    for (double value : values) {
        result = std::max(result, std::abs(value));
    }

    return result;
}

} // namespace

DynamicGibbsFlashCalculation::DynamicGibbsFlashCalculation(
    const DynamicPengRobinsonEOS& eos
    )
    : eos_(eos)
{
}

DynamicComposition DynamicGibbsFlashCalculation::normalizeComposition(
    const DynamicComposition& composition
    ) const {
    if (composition.size() != eos_.componentCount()) {
        throw std::runtime_error(
            "DynamicGibbsFlashCalculation: composition size mismatch"
            );
    }

    double sum = 0.0;

    for (double value : composition) {
        if (value < 0.0) {
            throw std::runtime_error(
                "DynamicGibbsFlashCalculation: negative composition value"
                );
        }

        sum += value;
    }

    if (sum <= 0.0) {
        throw std::runtime_error(
            "DynamicGibbsFlashCalculation: composition sum must be positive"
            );
    }

    DynamicComposition result = composition;

    for (double& value : result) {
        value /= sum;
    }

    return result;
}

DynamicComposition DynamicGibbsFlashCalculation::computeFugacityResiduals(
    double pressurePa,
    const DynamicComposition& xLiquid,
    const DynamicComposition& yVapor,
    const DynamicComposition& phiLiquid,
    const DynamicComposition& phiVapor
    ) const {
    const std::size_t n = eos_.componentCount();

    DynamicComposition residuals(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        const double fLiquid =
            safePositive(xLiquid[i] * phiLiquid[i] * pressurePa);

        const double fVapor =
            safePositive(yVapor[i] * phiVapor[i] * pressurePa);

        residuals[i] = std::log(fLiquid) - std::log(fVapor);
    }

    return residuals;
}

double DynamicGibbsFlashCalculation::computeGibbsDeviationDimensionless(
    const DynamicComposition& overallComposition,
    const DynamicComposition& liquidSplit,
    const DynamicGibbsFlashResult& state
    ) const {
    const std::size_t n = eos_.componentCount();

    double g = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        const double zi = overallComposition[i];

        if (zi <= 0.0) {
            continue;
        }

        const double nLiquid = zi * liquidSplit[i];
        const double nVapor = zi * (1.0 - liquidSplit[i]);

        const double liquidArgument =
            safePositive(state.xLiquid[i] * state.phiLiquid[i] / zi);

        const double vaporArgument =
            safePositive(state.yVapor[i] * state.phiVapor[i] / zi);

        g += nLiquid * std::log(liquidArgument);
        g += nVapor * std::log(vaporArgument);
    }

    return g;
}

DynamicGibbsFlashResult DynamicGibbsFlashCalculation::evaluateState(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& overallComposition,
    const DynamicComposition& liquidSplit
    ) const {
    const std::size_t n = eos_.componentCount();

    DynamicGibbsFlashResult result;

    result.liquidSplit = liquidSplit;

    result.xLiquid.assign(n, 0.0);
    result.yVapor.assign(n, 0.0);

    double betaLiquid = 0.0;
    double betaVapor = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        const double zi = overallComposition[i];

        const double nLiquid = zi * liquidSplit[i];
        const double nVapor = zi * (1.0 - liquidSplit[i]);

        betaLiquid += nLiquid;
        betaVapor += nVapor;
    }

    result.betaLiquid = betaLiquid;
    result.betaVapor = betaVapor;

    if (betaLiquid <= 0.0 || betaVapor <= 0.0) {
        result.singlePhase = true;
        result.status = "Degenerate phase split";
        return result;
    }

    for (std::size_t i = 0; i < n; ++i) {
        const double zi = overallComposition[i];

        result.xLiquid[i] =
            zi * liquidSplit[i] / betaLiquid;

        result.yVapor[i] =
            zi * (1.0 - liquidSplit[i]) / betaVapor;
    }

    result.zLiquid =
        eos_.computeZFactors(
            pressurePa,
            temperatureK,
            result.xLiquid
            );

    result.zVapor =
        eos_.computeZFactors(
            pressurePa,
            temperatureK,
            result.yVapor
            );

    result.phiLiquid =
        eos_.computeFugacityCoefficients(
            pressurePa,
            temperatureK,
            result.xLiquid,
            result.zLiquid.liquid
            );

    result.phiVapor =
        eos_.computeFugacityCoefficients(
            pressurePa,
            temperatureK,
            result.yVapor,
            result.zVapor.vapor
            );

    const DynamicComposition residuals =
        computeFugacityResiduals(
            pressurePa,
            result.xLiquid,
            result.yVapor,
            result.phiLiquid,
            result.phiVapor
            );

    result.maxFugacityResidual = maxAbsValue(residuals);

    result.gibbsDeviationDimensionless =
        computeGibbsDeviationDimensionless(
            overallComposition,
            liquidSplit,
            result
            );

    return result;
}

DynamicGibbsFlashResult DynamicGibbsFlashCalculation::calculate(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& overallComposition,
    const DynamicGibbsFlashOptions& options
    ) const {
    if (pressurePa <= 0.0 || temperatureK <= 0.0) {
        throw std::runtime_error(
            "DynamicGibbsFlashCalculation: pressure and temperature must be positive"
            );
    }

    const DynamicComposition z =
        normalizeComposition(overallComposition);

    const std::size_t n = eos_.componentCount();

    const DynamicZFactors feedZ =
        eos_.computeZFactors(
            pressurePa,
            temperatureK,
            z
            );

    if (feedZ.roots.size() < 2) {
        DynamicGibbsFlashResult result;

        result.singlePhase = true;
        result.converged = true;
        result.status = "Single phase by PR EOS roots";

        result.iterations = 0;

        result.betaLiquid = 0.0;
        result.betaVapor = 1.0;

        result.liquidSplit.assign(n, 0.0);

        result.xLiquid = z;
        result.yVapor = z;

        result.zLiquid = feedZ;
        result.zVapor = feedZ;

        result.phiLiquid =
            eos_.computeFugacityCoefficients(
                pressurePa,
                temperatureK,
                z,
                feedZ.vapor
                );

        result.phiVapor = result.phiLiquid;

        result.maxFugacityResidual = 0.0;
        result.gibbsDeviationDimensionless = 0.0;

        return result;
    }

    DynamicComposition split(
        n,
        clampValue(
            options.initialLiquidSplit,
            options.minSplit,
            1.0 - options.minSplit
            )
        );

    DynamicGibbsFlashResult state;

    for (int iteration = 0; iteration < options.maxIterations; ++iteration) {
        state =
            evaluateState(
                pressurePa,
                temperatureK,
                z,
                split
                );

        state.iterations = iteration;

        if (state.betaLiquid < options.minPhaseFraction) {
            state.singlePhase = true;
            state.converged = true;
            state.status = "Converged to vapor phase";
            return state;
        }

        if (state.betaVapor < options.minPhaseFraction) {
            state.singlePhase = true;
            state.converged = true;
            state.status = "Converged to liquid phase";
            return state;
        }

        if (state.maxFugacityResidual < options.tolerance) {
            state.converged = true;
            state.singlePhase = false;
            state.status = "Two-phase Gibbs ODE flash converged";
            return state;
        }

        const DynamicComposition residuals =
            computeFugacityResiduals(
                pressurePa,
                state.xLiquid,
                state.yVapor,
                state.phiLiquid,
                state.phiVapor
                );

        double probabilityFactor = 1.0;

        if (options.useGibbsProbabilityFactor) {
            probabilityFactor =
                std::exp(
                    -clampValue(
                        state.gibbsDeviationDimensionless,
                        -50.0,
                        50.0
                        )
                    );
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (z[i] <= 0.0) {
                split[i] = options.minSplit;
                continue;
            }

            const double boundFactor =
                split[i] * (1.0 - split[i]);

            double ds =
                -options.pseudoTimeStep
                * options.mobility
                * probabilityFactor
                * boundFactor
                * residuals[i];

            ds =
                clampValue(
                    ds,
                    -options.maxSplitStep,
                    options.maxSplitStep
                    );

            split[i] =
                clampValue(
                    split[i] + ds,
                    options.minSplit,
                    1.0 - options.minSplit
                    );
        }
    }

    state =
        evaluateState(
            pressurePa,
            temperatureK,
            z,
            split
            );

    state.converged = false;
    state.status = "Gibbs ODE flash did not converge";
    state.iterations = options.maxIterations;

    return state;
}

DynamicGibbsFlashResult
DynamicGibbsFlashCalculation::evaluateEquilibriumResidual(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& overallComposition,
    const DynamicComposition& liquidSplit
    ) const {
    const DynamicComposition z =
        normalizeComposition(overallComposition);

    return evaluateState(
        pressurePa,
        temperatureK,
        z,
        liquidSplit
        );
}

DynamicComposition DynamicGibbsFlashCalculation::fugacityResiduals(
    double pressurePa,
    const DynamicComposition& xLiquid,
    const DynamicComposition& yVapor,
    const DynamicComposition& phiLiquid,
    const DynamicComposition& phiVapor
    ) const {
    return computeFugacityResiduals(
        pressurePa,
        xLiquid,
        yVapor,
        phiLiquid,
        phiVapor
        );
}
