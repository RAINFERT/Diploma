#include "GibbsFlashCalculation.h"

#include "../numerics/NewtonSolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>


namespace
{
constexpr double Tiny = 1.0e-30;
constexpr double ActiveTolerance = 1.0e-14;
constexpr double PhaseFractionMin = 1.0e-10;
constexpr double MoleFloorRelative = 1.0e-12;
constexpr double MinPseudoDriftScale = 1.0e-8;

constexpr double BarrierOmega = 1.0e5;
constexpr double BarrierWidthFactor = 0.1;
constexpr double Pi = 3.14159265358979323846;

double safeLog(double value)
{
    return std::log(std::max(value, Tiny));
}


double phaseMoleFloor(double zi)
{
    return std::max(Tiny, MoleFloorRelative * std::max(zi, Tiny));
}

double barrierWidth(double zi)
{
    // Ширина сглаживания должна быть намного меньше допустимого
    // количества вещества компонента в фазе.
    return std::max(1.0e-18, BarrierWidthFactor * phaseMoleFloor(zi));
}

double smoothHeaviside(double x, double zi)
{
    const double width = barrierWidth(zi);
    return 0.5 * (1.0 + std::erf(x / width));
}

double smoothDelta(double x, double zi)
{
    const double width = barrierWidth(zi);
    const double q = x / width;

    return std::exp(-0.5 * q * q)
           / (std::sqrt(2.0 * Pi) * width);
}

double clampExpArgument(double value)
{
    return std::clamp(value, -50.0, 50.0);
}

double sumComposition(const Composition& composition)
{
    double sum = 0.0;

    for (double value : composition)
    {
        sum += value;
    }

    return sum;
}
}


GibbsFlashCalculation::GibbsFlashCalculation(
    const MaterialList& materials,
    const PengRobinsonEOS& eos
    )
    : materials_(materials),
    eos_(eos)
{
}


bool GibbsFlashCalculation::isActive(
    double zi
    ) const
{
    return zi > ActiveTolerance;
}


Composition GibbsFlashCalculation::normalizeComposition(
    const Composition& composition
    ) const
{
    const double sum = sumComposition(composition);

    if (sum <= 0.0)
    {
        throw std::runtime_error("Overall composition sum must be positive");
    }

    Composition result{};

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (composition[i] < 0.0)
        {
            throw std::runtime_error("Composition contains negative value");
        }

        result[i] = composition[i] / sum;
    }

    return result;
}


Composition GibbsFlashCalculation::wilsonKValues(
    double pressurePa,
    double temperatureK
    ) const
{
    Composition k{};

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const Material& material = materials_[i];

        k[i] =
            material.criticalPressurePa / pressurePa
            *
            std::exp(
                5.373
                * (1.0 + material.acentricFactor)
                * (1.0 - material.criticalTemperatureK / temperatureK)
                );

        k[i] = std::clamp(k[i], 1.0e-12, 1.0e12);
    }

    return k;
}


double GibbsFlashCalculation::rachfordRiceValue(
    double beta,
    const Composition& z,
    const Composition& k
    ) const
{
    double value = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (!isActive(z[i]))
        {
            continue;
        }

        const double denominator =
            1.0 + beta * (k[i] - 1.0);

        value += z[i] * (k[i] - 1.0) / denominator;
    }

    return value;
}


double GibbsFlashCalculation::initialBetaEstimate(
    const Composition& z,
    const Composition& k
    ) const
{
    const double f0 = rachfordRiceValue(0.0, z, k);
    const double f1 = rachfordRiceValue(1.0, z, k);

    if (f0 <= 0.0)
    {
        return 1.0e-4;
    }

    if (f1 >= 0.0)
    {
        return 1.0 - 1.0e-4;
    }

    double left = 0.0;
    double right = 1.0;

    for (int iter = 0; iter < 100; ++iter)
    {
        const double mid = 0.5 * (left + right);
        const double fmid = rachfordRiceValue(mid, z, k);

        if (fmid > 0.0)
        {
            left = mid;
        }
        else
        {
            right = mid;
        }
    }

    return std::clamp(0.5 * (left + right), 1.0e-4, 1.0 - 1.0e-4);
}

GibbsFlashCalculation::GibbsState GibbsFlashCalculation::initialStateFromBeta(
    double pressurePa,
    double temperatureK,
    const Composition& z,
    const Composition& k,
    double betaVapor
    ) const
{
    const double betaV =
        std::clamp(betaVapor, 1.0e-10, 1.0 - 1.0e-10);

    Composition x{};
    Composition y{};

    double sumX = 0.0;
    double sumY = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (!isActive(z[i]))
        {
            x[i] = 0.0;
            y[i] = 0.0;
            continue;
        }

        const double denominator =
            1.0 + betaV * (k[i] - 1.0);

        if (denominator <= 0.0)
        {
            throw std::runtime_error(
                "Invalid denominator in initialStateFromBeta"
                );
        }

        x[i] = z[i] / denominator;
        y[i] = k[i] * x[i];

        sumX += x[i];
        sumY += y[i];
    }

    if (sumX <= 0.0 || sumY <= 0.0)
    {
        throw std::runtime_error(
            "Invalid initial phase compositions"
            );
    }

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (!isActive(z[i]))
        {
            continue;
        }

        x[i] /= sumX;
        y[i] /= sumY;
    }

    double vLiquid =
        R * temperatureK / pressurePa * 0.02;

    double vVapor =
        R * temperatureK / pressurePa;

    try
    {
        const ZFactors zL =
            eos_.computeZFactors(
                pressurePa,
                temperatureK,
                x
                );

        const ZFactors zV =
            eos_.computeZFactors(
                pressurePa,
                temperatureK,
                y
                );

        vLiquid =
            zL.liquid * R * temperatureK / pressurePa;

        vVapor =
            zV.vapor * R * temperatureK / pressurePa;
    }
    catch (...)
    {
        // Оставляем безопасную оценку.
    }

    const double idealGasVolume =
        R * temperatureK / pressurePa;

    if (vLiquid <= 0.0 || !std::isfinite(vLiquid))
    {
        vLiquid = 0.02 * idealGasVolume;
    }

    if (vVapor <= vLiquid || !std::isfinite(vVapor))
    {
        vVapor = idealGasVolume;
    }

    GibbsState state;
    state.y.resize(2 + ComponentCount);

    state.y[0] = std::log(vLiquid);
    state.y[1] = std::log(vVapor);

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (!isActive(z[i]))
        {
            state.y[2 + i] = std::log(Tiny);
            continue;
        }

        double nLiquid =
            (1.0 - betaV) * x[i];

        const double floor = phaseMoleFloor(z[i]);
        const double upper = std::max(z[i] - floor, floor);

        nLiquid =
            std::clamp(nLiquid, floor, upper);

        state.y[2 + i] =
            std::log(nLiquid);
    }

    return state;
}

GibbsFlashCalculation::GibbsState GibbsFlashCalculation::initialState(
    double pressurePa,
    double temperatureK,
    const Composition& z
    ) const
{
    const Composition k =
        wilsonKValues(
            pressurePa,
            temperatureK
            );

    std::vector<double> betaCandidates;

    betaCandidates.push_back(
        initialBetaEstimate(z, k)
        );

    // Особенно важны малые beta для случаев вроде:
    // почти чистая вода + следы летучего компонента.
    betaCandidates.push_back(1.0e-10);
    betaCandidates.push_back(1.0e-8);
    betaCandidates.push_back(1.0e-6);
    betaCandidates.push_back(1.0e-5);
    betaCandidates.push_back(1.0e-4);
    betaCandidates.push_back(3.0e-4);
    betaCandidates.push_back(1.0e-3);
    betaCandidates.push_back(3.0e-3);
    betaCandidates.push_back(1.0e-2);
    betaCandidates.push_back(3.0e-2);
    betaCandidates.push_back(1.0e-1);
    betaCandidates.push_back(3.0e-1);
    betaCandidates.push_back(5.0e-1);
    betaCandidates.push_back(7.0e-1);
    betaCandidates.push_back(9.0e-1);

    bool found = false;
    double bestScore =
        std::numeric_limits<double>::infinity();

    GibbsState bestState;

    for (double beta : betaCandidates)
    {
        beta =
            std::clamp(
                beta,
                1.0e-10,
                1.0 - 1.0e-10
                );

        try
        {
            GibbsState candidate =
                initialStateFromBeta(
                    pressurePa,
                    temperatureK,
                    z,
                    k,
                    beta
                    );

            const GibbsEvaluation e =
                evaluate(
                    candidate,
                    pressurePa,
                    temperatureK,
                    z
                    );

            if (!std::isfinite(e.gibbsDimensionless))
            {
                continue;
            }

            // Это не критерий равновесия.
            // Это только выбор лучшего начального приближения.
            const double score =
                e.gibbsDimensionless
                + 1.0e-3 * e.maxFugacityResidual
                + 1.0e-3 * e.maxPressureResidual;

            if (score < bestScore)
            {
                bestScore = score;
                bestState = candidate;
                found = true;
            }
        }
        catch (...)
        {
            // Неподходящий старт просто пропускаем.
        }
    }

    if (!found)
    {
        return initialStateFromBeta(
            pressurePa,
            temperatureK,
            z,
            k,
            0.5
            );
    }

    return bestState;
}

GibbsFlashCalculation::GibbsEvaluation GibbsFlashCalculation::evaluate(
    const GibbsState& state,
    double pressurePa,
    double temperatureK,
    const Composition& z
    ) const
{
    GibbsEvaluation e;

    e.vLiquid = std::exp(state.y[0]);
    e.vVapor = std::exp(state.y[1]);

    if (e.vLiquid <= 0.0 || e.vVapor <= 0.0)
    {
        throw std::runtime_error("Invalid molar volume in Gibbs state");
    }

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (!isActive(z[i]))
        {
            e.nLiquid[i] = 0.0;
            e.nVapor[i] = 0.0;
            continue;
        }

        const double nL = std::exp(state.y[2 + i]);
        const double floor = phaseMoleFloor(z[i]);

        if (nL <= floor || nL >= z[i] - floor)
        {
            throw std::runtime_error("Invalid phase mole amount in Gibbs state");
        }

        e.nLiquid[i] = nL;
        e.nVapor[i] = z[i] - nL;

        e.betaLiquid += e.nLiquid[i];
        e.betaVapor += e.nVapor[i];
    }

    if (e.betaLiquid <= PhaseFractionMin || e.betaVapor <= PhaseFractionMin)
    {
        throw std::runtime_error("Degenerate phase split in Gibbs state");
    }

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (!isActive(z[i]))
        {
            e.xLiquid[i] = 0.0;
            e.yVapor[i] = 0.0;
            continue;
        }

        e.xLiquid[i] = e.nLiquid[i] / e.betaLiquid;
        e.yVapor[i] = e.nVapor[i] / e.betaVapor;
    }

    e.zLiquid = pressurePa * e.vLiquid / (R * temperatureK);
    e.zVapor = pressurePa * e.vVapor / (R * temperatureK);

    e.lnPhiLiquid =
        eos_.computeLogFugacityCoefficients(
            pressurePa,
            temperatureK,
            e.xLiquid,
            e.zLiquid
            );

    e.lnPhiVapor =
        eos_.computeLogFugacityCoefficients(
            pressurePa,
            temperatureK,
            e.yVapor,
            e.zVapor
            );

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        e.phiLiquid[i] = std::exp(e.lnPhiLiquid[i]);
        e.phiVapor[i] = std::exp(e.lnPhiVapor[i]);
    }

    e.pressureLiquidPa =
        eos_.pressureFromMolarVolume(
            temperatureK,
            e.vLiquid,
            e.xLiquid
            );

    e.pressureVaporPa =
        eos_.pressureFromMolarVolume(
            temperatureK,
            e.vVapor,
            e.yVapor
            );

    double g = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (!isActive(z[i]))
        {
            continue;
        }

        const double liquidTerm =
            e.nLiquid[i]
            * (
                safeLog(e.nLiquid[i])
                + e.lnPhiLiquid[i]
                - safeLog(e.betaLiquid)
                );

        const double vaporTerm =
            e.nVapor[i]
            * (
                safeLog(e.nVapor[i])
                + e.lnPhiVapor[i]
                - safeLog(e.betaVapor)
                );

        const double barrier =
            BarrierOmega
            * smoothHeaviside(e.nLiquid[i] - z[i], z[i]);

        g += liquidTerm + vaporTerm + barrier;
    }

    e.gibbsDimensionless = g;

    e.maxFugacityResidual = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (!isActive(z[i]))
        {
            continue;
        }

        const double residual =
            safeLog(e.xLiquid[i])
            + e.lnPhiLiquid[i]
            - safeLog(e.yVapor[i])
            - e.lnPhiVapor[i];

        e.maxFugacityResidual =
            std::max(e.maxFugacityResidual, std::abs(residual));
    }

    const double liquidPressureResidual =
        std::abs(e.pressureLiquidPa - pressurePa) / pressurePa;

    const double vaporPressureResidual =
        std::abs(e.pressureVaporPa - pressurePa) / pressurePa;

    e.maxPressureResidual =
        std::max(liquidPressureResidual, vaporPressureResidual);

    return e;
}


std::vector<double> GibbsFlashCalculation::rhs(
    const GibbsState& state,
    double pressurePa,
    double temperatureK,
    const Composition& z
    ) const
{
    const GibbsEvaluation e =
        evaluate(state, pressurePa, temperatureK, z);

    std::vector<double> dy(2 + ComponentCount, 0.0);

    const double vcMix =
        eos_.pseudoCriticalMolarVolumeM3PerMol(z);

    const double scale =
        std::max(
            MinPseudoDriftScale,
            std::exp(clampExpArgument(-e.gibbsDimensionless))
            );

    dy[0] =
        scale
        * vcMix
        * (e.pressureLiquidPa - pressurePa)
        / (R * temperatureK)
        / e.vLiquid;

    dy[1] =
        scale
        * vcMix
        * (e.pressureVaporPa - pressurePa)
        / (R * temperatureK)
        / e.vVapor;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (!isActive(z[i]))
        {
            dy[2 + i] = 0.0;
            continue;
        }

        const double liquidTerm =
            safeLog(e.nLiquid[i])
            + e.lnPhiLiquid[i]
            - safeLog(e.betaLiquid);

        const double vaporTerm =
            safeLog(e.nVapor[i])
            + e.lnPhiVapor[i]
            - safeLog(e.betaVapor);

        const double barrierDerivative =
            BarrierOmega
            * e.nLiquid[i]
            * smoothDelta(e.nLiquid[i] - z[i], z[i]);

        dy[2 + i] =
            scale
            * (
                vaporTerm
                - liquidTerm
                - barrierDerivative
                );
    }

    return dy;
}

GibbsFlashCalculation::RadauStepResult GibbsFlashCalculation::radauIIA3Step(
    const GibbsState& state,
    double pseudoTimeStep,
    double pressurePa,
    double temperatureK,
    const Composition& z
    ) const
{
    using numerics::Vector;

    const std::size_t n = state.y.size();
    const std::size_t stages = 3;

    const double sqrt6 = std::sqrt(6.0);

    const double A[3][3] =
        {
            {
                (88.0 - 7.0 * sqrt6) / 360.0,
                (296.0 - 169.0 * sqrt6) / 1800.0,
                (-2.0 + 3.0 * sqrt6) / 225.0
            },
            {
                (296.0 + 169.0 * sqrt6) / 1800.0,
                (88.0 + 7.0 * sqrt6) / 360.0,
                (-2.0 - 3.0 * sqrt6) / 225.0
            },
            {
                (16.0 - sqrt6) / 36.0,
                (16.0 + sqrt6) / 36.0,
                1.0 / 9.0
            }
        };

    Vector initialGuess(stages * n, 0.0);

    for (std::size_t s = 0; s < stages; ++s)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            initialGuess[s * n + i] = state.y[i];
        }
    }

    numerics::NewtonOptions options;
    options.maxIterations = 40;
    options.residualTolerance = 1.0e-9;
    options.stepTolerance = 1.0e-10;
    options.finiteDifferenceStep = 1.0e-6;
    options.minDampingFactor = 1.0e-8;

    options.lowerBounds.resize(stages * n);
    options.upperBounds.resize(stages * n);

    const double idealGasVolume =
        R * temperatureK / pressurePa;

    const double minVolume =
        std::max(1.0e-10, 1.0e-6 * idealGasVolume);

    const double maxVolume =
        std::max(1.0e-4, 1.0e3 * idealGasVolume);

    for (std::size_t s = 0; s < stages; ++s)
    {
        options.lowerBounds[s * n + 0] = std::log(minVolume);
        options.upperBounds[s * n + 0] = std::log(maxVolume);

        options.lowerBounds[s * n + 1] = std::log(minVolume);
        options.upperBounds[s * n + 1] = std::log(maxVolume);

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            const std::size_t index = s * n + 2 + i;

            if (!isActive(z[i]))
            {
                const double fixed = std::log(Tiny);

                options.lowerBounds[index] = fixed;
                options.upperBounds[index] = fixed;

                continue;
            }

            const double floor = phaseMoleFloor(z[i]);
            const double upper = std::max(z[i] - floor, floor);

            options.lowerBounds[index] = std::log(floor);
            options.upperBounds[index] = std::log(upper);
        }
    }
    int rhsEvaluations = 0;
    auto residualFunction =
        [&](const Vector& unknown) -> Vector
    {
        std::array<GibbsState, 3> stageStates;
        std::array<Vector, 3> stageRhs;

        for (std::size_t s = 0; s < stages; ++s)
        {
            stageStates[s].y.resize(n);

            for (std::size_t i = 0; i < n; ++i)
            {
                stageStates[s].y[i] = unknown[s * n + i];
            }

            stageRhs[s] =
                rhs(
                    stageStates[s],
                    pressurePa,
                    temperatureK,
                    z
                    );

            ++rhsEvaluations;
        }

        Vector residual(stages * n, 0.0);

        for (std::size_t s = 0; s < stages; ++s)
        {
            for (std::size_t i = 0; i < n; ++i)
            {
                double collocationValue = state.y[i];

                for (std::size_t m = 0; m < stages; ++m)
                {
                    collocationValue +=
                        pseudoTimeStep * A[s][m] * stageRhs[m][i];
                }

                residual[s * n + i] =
                    stageStates[s].y[i] - collocationValue;
            }
        }

        return residual;
    };

    const numerics::NewtonResult result =
        numerics::NewtonSolver::solve(
            residualFunction,
            initialGuess,
            options
            );

    if (!result.converged)
    {
        throw std::runtime_error(
            "Gibbs flash Radau step failed: " + result.message
            );
    }

    RadauStepResult stepResult;

    stepResult.state.y.resize(n);

    for (std::size_t i = 0; i < n; ++i)
    {
        stepResult.state.y[i] = result.solution[2 * n + i];
    }

    stepResult.newtonIterations = result.iterations;
    stepResult.newtonResidualNorm = result.residualNorm;
    stepResult.newtonStepNorm = result.stepNorm;
    stepResult.rhsEvaluations = rhsEvaluations;

    return stepResult;
}

double GibbsFlashCalculation::singlePhaseGibbsDimensionless(
    double pressurePa,
    double temperatureK,
    const Composition& z,
    bool vaporPhase
    ) const
{
    const ZFactors zFactors =
        eos_.computeZFactors(
            pressurePa,
            temperatureK,
            z
            );

    const double zFactor =
        vaporPhase ? zFactors.vapor : zFactors.liquid;

    const Composition lnPhi =
        eos_.computeLogFugacityCoefficients(
            pressurePa,
            temperatureK,
            z,
            zFactor
            );

    double g = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (!isActive(z[i]))
        {
            continue;
        }

        g += z[i] * (safeLog(z[i]) + lnPhi[i]);
    }

    return g;
}


FlashStatus GibbsFlashCalculation::bestSinglePhase(
    double pressurePa,
    double temperatureK,
    const Composition& z,
    double& bestGibbs
    ) const
{
    const double gLiquid =
        singlePhaseGibbsDimensionless(
            pressurePa,
            temperatureK,
            z,
            false
            );

    const double gVapor =
        singlePhaseGibbsDimensionless(
            pressurePa,
            temperatureK,
            z,
            true
            );

    if (gLiquid <= gVapor)
    {
        bestGibbs = gLiquid;
        return FlashStatus::SingleLiquid;
    }

    bestGibbs = gVapor;
    return FlashStatus::SingleVapor;
}


bool GibbsFlashCalculation::isConverged(
    const GibbsEvaluation& evaluation,
    double tolerance
    ) const
{
    return
        evaluation.maxFugacityResidual < tolerance
        &&
        evaluation.maxPressureResidual < tolerance;
}

FlashResult GibbsFlashCalculation::calculate(
    double pressurePa,
    double temperatureK,
    const Composition& zOverall,
    double tolerance,
    int maxIterations
    ) const
{
    if (pressurePa <= 0.0)
    {
        throw std::invalid_argument("Pressure must be positive");
    }

    if (temperatureK <= 0.0)
    {
        throw std::invalid_argument("Temperature must be positive");
    }

    const Composition z =
        normalizeComposition(zOverall);

    FlashDiagnostics diagnostics;

    const RachfordRicePrecheckResult precheck =
        rachfordRicePrecheck(
            pressurePa,
            temperatureK,
            z
            );

    diagnostics.precheckBetaVapor =
        precheck.betaVapor;

    diagnostics.precheckF0 =
        precheck.f0;

    diagnostics.precheckF1 =
        precheck.f1;

    diagnostics.precheckIterations =
        precheck.iterations;

    diagnostics.gibbsSingleLiquid =
        singlePhaseGibbsDimensionless(
            pressurePa,
            temperatureK,
            z,
            false
            );

    diagnostics.gibbsSingleVapor =
        singlePhaseGibbsDimensionless(
            pressurePa,
            temperatureK,
            z,
            true
            );

    const FlashStatus bestSingle =
        diagnostics.gibbsSingleLiquid <= diagnostics.gibbsSingleVapor
            ? FlashStatus::SingleLiquid
            : FlashStatus::SingleVapor;

    const double bestSingleGibbs =
        std::min(
            diagnostics.gibbsSingleLiquid,
            diagnostics.gibbsSingleVapor
            );

    if (precheck.action == PrecheckAction::SingleLiquid)
    {
        FlashResult result =
            makeSinglePhaseResult(
                FlashStatus::SingleLiquid,
                pressurePa,
                temperatureK,
                z
                );

        result.method =
            FlashMethod::HybridSinglePhasePrecheck;

        result.diagnostics =
            diagnostics;

        result.iterations = 0;

        return result;
    }

    if (precheck.action == PrecheckAction::SingleVapor)
    {
        FlashResult result =
            makeSinglePhaseResult(
                FlashStatus::SingleVapor,
                pressurePa,
                temperatureK,
                z
                );

        result.method =
            FlashMethod::HybridSinglePhasePrecheck;

        result.diagnostics =
            diagnostics;

        result.iterations = 0;

        return result;
    }

    if (precheck.action == PrecheckAction::BoundaryRachfordRice)
    {
        return makeBoundaryRachfordRiceResult(
            pressurePa,
            temperatureK,
            z,
            tolerance,
            maxIterations,
            diagnostics
            );
    }

    GibbsState state =
        initialState(
            pressurePa,
            temperatureK,
            z
            );

    try
    {
        const GibbsEvaluation initialEvaluation =
            evaluate(
                state,
                pressurePa,
                temperatureK,
                z
                );

        diagnostics.initialBetaVapor =
            initialEvaluation.betaVapor;

        diagnostics.initialGibbsTwoPhase =
            initialEvaluation.gibbsDimensionless;

        diagnostics.initialMaxFugacityResidual =
            initialEvaluation.maxFugacityResidual;

        diagnostics.initialMaxPressureResidual =
            initialEvaluation.maxPressureResidual;
    }
    catch (const std::exception& ex)
    {
        diagnostics.lastRadauFailureMessage =
            std::string("Initial evaluation failed: ") + ex.what();
    }
    catch (...)
    {
        diagnostics.lastRadauFailureMessage =
            "Initial evaluation failed: unknown exception";
    }

    GibbsEvaluation lastEvaluation;

    double pseudoTimeStep = 1.0e-2;

    for (int iter = 0; iter < maxIterations; ++iter)
    {
        try
        {
            lastEvaluation =
                evaluate(
                    state,
                    pressurePa,
                    temperatureK,
                    z
                    );

            diagnostics.betaLiquidCandidate = lastEvaluation.betaLiquid;
            diagnostics.betaVaporCandidate = lastEvaluation.betaVapor;

            diagnostics.pseudoSteps = iter + 1;
            diagnostics.gibbsTwoPhase = lastEvaluation.gibbsDimensionless;
            diagnostics.maxFugacityResidual = lastEvaluation.maxFugacityResidual;
            diagnostics.maxPressureResidual = lastEvaluation.maxPressureResidual;
            diagnostics.finalPseudoTimeStep = pseudoTimeStep;

            if (isConverged(lastEvaluation, tolerance))
            {
                const bool goodPhaseFractions =
                    lastEvaluation.betaVapor > PhaseFractionMin
                    &&
                    lastEvaluation.betaLiquid > PhaseFractionMin;

                const bool twoPhaseBetter =
                    lastEvaluation.gibbsDimensionless
                    <
                    bestSingleGibbs - 1.0e-10;

                if (goodPhaseFractions && twoPhaseBetter)
                {
                    FlashResult result =
                        makeTwoPhaseResult(
                            lastEvaluation,
                            iter
                            );

                    diagnostics.gibbsTwoPhase =
                        lastEvaluation.gibbsDimensionless;

                    diagnostics.maxFugacityResidual =
                        lastEvaluation.maxFugacityResidual;

                    diagnostics.maxPressureResidual =
                        lastEvaluation.maxPressureResidual;

                    diagnostics.finalPseudoTimeStep =
                        pseudoTimeStep;

                    result.diagnostics = diagnostics;

                    fillDensities(
                        result,
                        pressurePa,
                        temperatureK
                        );

                    return result;
                }

                FlashResult result =
                    makeSinglePhaseResult(
                        bestSingle,
                        pressurePa,
                        temperatureK,
                        z
                        );

                result.diagnostics = diagnostics;

                return result;
            }

            bool stepAccepted = false;

            for (int attempt = 0; attempt < 8; ++attempt)
            {
                try
                {
                    RadauStepResult radauStep =
                        radauIIA3Step(
                            state,
                            pseudoTimeStep,
                            pressurePa,
                            temperatureK,
                            z
                            );

                    GibbsState candidate = radauStep.state;

                    const GibbsEvaluation candidateEvaluation =
                        evaluate(
                            candidate,
                            pressurePa,
                            temperatureK,
                            z
                            );

                    diagnostics.betaLiquidCandidate = lastEvaluation.betaLiquid;
                    diagnostics.betaVaporCandidate = lastEvaluation.betaVapor;

                    if (
                        std::isfinite(candidateEvaluation.gibbsDimensionless)
                        &&
                        candidateEvaluation.maxPressureResidual < 1.0e6
                        )
                    {
                        diagnostics.acceptedRadauSteps += 1;

                        diagnostics.totalRadauNewtonIterations +=
                            radauStep.newtonIterations;

                        diagnostics.lastRadauNewtonIterations =
                            radauStep.newtonIterations;

                        diagnostics.lastRadauResidualNorm =
                            radauStep.newtonResidualNorm;

                        diagnostics.lastRadauStepNorm =
                            radauStep.newtonStepNorm;

                        diagnostics.rhsEvaluations +=
                            radauStep.rhsEvaluations;

                        state = candidate;

                        pseudoTimeStep =
                            std::min(1.0, pseudoTimeStep * 1.2);

                        stepAccepted = true;
                        break;
                    }

                    diagnostics.rejectedRadauSteps += 1;
                    pseudoTimeStep *= 0.5;
                }
                catch (const std::exception& ex)
                {
                    pseudoTimeStep *= 0.5;
                    diagnostics.rejectedRadauSteps += 1;
                    diagnostics.lastRadauFailureMessage = ex.what();
                }
                catch (...)
                {
                    pseudoTimeStep *= 0.5;
                    diagnostics.rejectedRadauSteps += 1;
                    diagnostics.lastRadauFailureMessage =
                        "Unknown Radau failure";
                }
            }

            if (!stepAccepted)
            {
                break;
            }
        }
        catch (...)
        {
            break;
        }
    }

    // Если двухфазный Radau не сошелся — возвращаем лучший однофазный вариант.
    FlashResult result =
        makeSinglePhaseResult(
            bestSingle,
            pressurePa,
            temperatureK,
            z
            );

    diagnostics.finalPseudoTimeStep = pseudoTimeStep;

    result.diagnostics = diagnostics;

    // Для fallback-случая показываем не только принятые шаги,
    // а общее число попыток Radau.
    result.iterations =
        diagnostics.acceptedRadauSteps
        + diagnostics.rejectedRadauSteps;

    return result;
}

FlashResult GibbsFlashCalculation::makeTwoPhaseResult(
    const GibbsEvaluation& evaluation,
    int iterations
    ) const
{
    FlashResult result;

    result.method = FlashMethod::GibbsRadau;
    result.status = FlashStatus::Converged;
    result.iterations = iterations;

    result.beta = evaluation.betaVapor;

    result.xLiquid = evaluation.xLiquid;
    result.yVapor = evaluation.yVapor;

    result.zLiquid = evaluation.zLiquid;
    result.zVapor = evaluation.zVapor;

    result.phiLiquid = evaluation.phiLiquid;
    result.phiVapor = evaluation.phiVapor;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (result.xLiquid[i] > 0.0)
        {
            result.kValues[i] =
                result.yVapor[i] / result.xLiquid[i];
        }
        else
        {
            result.kValues[i] = 0.0;
        }
    }

    return result;
}


FlashResult GibbsFlashCalculation::makeSinglePhaseResult(
    FlashStatus status,
    double pressurePa,
    double temperatureK,
    const Composition& z
    ) const
{
    FlashResult result;

    result.method = FlashMethod::GibbsRadau;
    result.status = status;
    result.iterations = 0;

    result.xLiquid = z;
    result.yVapor = z;

    const ZFactors zFactors =
        eos_.computeZFactors(
            pressurePa,
            temperatureK,
            z
            );

    if (status == FlashStatus::SingleVapor)
    {
        result.beta = 1.0;
        result.zVapor = zFactors.vapor;
        result.zLiquid = zFactors.vapor;

        result.phiVapor =
            eos_.computeFugacityCoefficients(
                pressurePa,
                temperatureK,
                z,
                result.zVapor
                );

        result.phiLiquid = result.phiVapor;
    }
    else
    {
        result.beta = 0.0;
        result.zLiquid = zFactors.liquid;
        result.zVapor = zFactors.liquid;

        result.phiLiquid =
            eos_.computeFugacityCoefficients(
                pressurePa,
                temperatureK,
                z,
                result.zLiquid
                );

        result.phiVapor = result.phiLiquid;
    }

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        result.kValues[i] = 1.0;
    }

    fillDensities(result, pressurePa, temperatureK);

    return result;
}


double GibbsFlashCalculation::averageMolarMassKgPerKmol(
    const Composition& composition
    ) const
{
    double value = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        value +=
            composition[i]
            * materials_[i].molarMassKgPerKmol;
    }

    return value;
}


void GibbsFlashCalculation::fillDensities(
    FlashResult& result,
    double pressurePa,
    double temperatureK
    ) const
{
    result.molarMassLiquidKgPerKmol =
        averageMolarMassKgPerKmol(result.xLiquid);

    result.molarMassVaporKgPerKmol =
        averageMolarMassKgPerKmol(result.yVapor);

    result.molarDensityLiquidKmolPerM3 =
        pressurePa
        / (result.zLiquid * R * temperatureK)
        / 1000.0;

    result.molarDensityVaporKmolPerM3 =
        pressurePa
        / (result.zVapor * R * temperatureK)
        / 1000.0;

    result.massDensityLiquidKgPerM3 =
        result.molarDensityLiquidKmolPerM3
        * result.molarMassLiquidKgPerKmol;

    result.massDensityVaporKgPerM3 =
        result.molarDensityVaporKmolPerM3
        * result.molarMassVaporKgPerKmol;

    const double beta = result.beta;

    result.molarDensityMixtureKmolPerM3 =
        1.0
        /
        (
            (1.0 - beta)
                / std::max(result.molarDensityLiquidKmolPerM3, 1.0e-30)
            +
            beta
                / std::max(result.molarDensityVaporKmolPerM3, 1.0e-30)
            );

    const double mixtureMolarMass =
        (1.0 - beta) * result.molarMassLiquidKgPerKmol
        +
        beta * result.molarMassVaporKgPerKmol;

    result.massDensityMixtureKgPerM3 =
        result.molarDensityMixtureKmolPerM3
        * mixtureMolarMass;
}

GibbsFlashCalculation::RachfordRicePrecheckResult
GibbsFlashCalculation::rachfordRicePrecheck(
    double pressurePa,
    double temperatureK,
    const Composition& z
    ) const
{
    constexpr double BoundaryBetaMin = 0.01;
    constexpr double BoundaryBetaMax = 0.99;

    RachfordRicePrecheckResult result;

    result.kValues =
        wilsonKValues(
            pressurePa,
            temperatureK
            );

    result.f0 =
        rachfordRiceValue(
            0.0,
            z,
            result.kValues
            );

    result.f1 =
        rachfordRiceValue(
            1.0,
            z,
            result.kValues
            );

    // Условие существования корня Рачфорда–Райса:
    //
    // F(0) > 0 и F(1) < 0.
    //
    // Если F(0) <= 0, смесь находится в жидкой области.
    // Если F(1) >= 0, смесь находится в паровой области.
    if (result.f0 <= 0.0)
    {
        result.action = PrecheckAction::SingleLiquid;
        result.betaVapor = 0.0;
        return result;
    }

    if (result.f1 >= 0.0)
    {
        result.action = PrecheckAction::SingleVapor;
        result.betaVapor = 1.0;
        return result;
    }

    result.betaVapor =
        solveRachfordRiceBeta(
            z,
            result.kValues,
            result.iterations
            );

    if (
        result.betaVapor < BoundaryBetaMin
        ||
        result.betaVapor > BoundaryBetaMax
        )
    {
        result.action = PrecheckAction::BoundaryRachfordRice;
        return result;
    }

    result.action = PrecheckAction::RunGibbsRadau;
    return result;
}

double GibbsFlashCalculation::solveRachfordRiceBeta(
    const Composition& z,
    const Composition& k,
    int& iterations
    ) const
{
    double left = 0.0;
    double right = 1.0;

    iterations = 0;

    double fLeft =
        rachfordRiceValue(
            left,
            z,
            k
            );

    double fRight =
        rachfordRiceValue(
            right,
            z,
            k
            );

    if (fLeft <= 0.0)
    {
        return 0.0;
    }

    if (fRight >= 0.0)
    {
        return 1.0;
    }

    for (int iter = 0; iter < 100; ++iter)
    {
        iterations = iter + 1;

        const double mid =
            0.5 * (left + right);

        const double fMid =
            rachfordRiceValue(
                mid,
                z,
                k
                );

        if (std::abs(fMid) < 1.0e-12)
        {
            return mid;
        }

        if (fMid > 0.0)
        {
            left = mid;
        }
        else
        {
            right = mid;
        }

        if (std::abs(right - left) < 1.0e-12)
        {
            break;
        }
    }

    return 0.5 * (left + right);
}

FlashResult GibbsFlashCalculation::makeBoundaryRachfordRiceResult(
    double pressurePa,
    double temperatureK,
    const Composition& z,
    double tolerance,
    int maxIterations,
    const FlashDiagnostics& diagnostics
    ) const
{
    FlashCalculation flashCalculation(
        materials_,
        eos_
        );

    FlashResult result =
        flashCalculation.calculate(
            pressurePa,
            temperatureK,
            z,
            tolerance,
            std::max(maxIterations, 80)
            );

    result.method =
        FlashMethod::HybridBoundaryRachfordRice;

    result.diagnostics =
        diagnostics;

    return result;
}
