#pragma once

#include "FlashCalculation.h"
#include "Material.h"
#include "PengRobinsonEOS.h"

#include <array>
#include <vector>


class GibbsFlashCalculation
{
public:
    GibbsFlashCalculation(
        const MaterialList& materials,
        const PengRobinsonEOS& eos
    );

    FlashResult calculate(
        double pressurePa,
        double temperatureK,
        const Composition& zOverall,
        double tolerance = 1.0e-8,
        int maxIterations = 300
    ) const;

private:
    static constexpr double R = 8.314462618;

    struct GibbsState
    {
        // y[0] = ln(V_liq), V_liq in m3/mol
        // y[1] = ln(V_vap), V_vap in m3/mol
        // y[2+i] = ln(n_liq_i)
        //
        // Общее количество смеси принято равным 1 моль:
        // z_i = n_liq_i + n_vap_i.
        std::vector<double> y;
    };

    struct RadauStepResult
    {
        GibbsState state;

        int newtonIterations = 0;
        double newtonResidualNorm = 0.0;
        double newtonStepNorm = 0.0;

        int rhsEvaluations = 0;
    };

    struct GibbsEvaluation
    {
        double vLiquid = 0.0;
        double vVapor = 0.0;

        Composition nLiquid{};
        Composition nVapor{};

        double betaLiquid = 0.0;
        double betaVapor = 0.0;

        Composition xLiquid{};
        Composition yVapor{};

        double zLiquid = 0.0;
        double zVapor = 0.0;

        Composition lnPhiLiquid{};
        Composition lnPhiVapor{};

        Composition phiLiquid{};
        Composition phiVapor{};

        double pressureLiquidPa = 0.0;
        double pressureVaporPa = 0.0;

        double gibbsDimensionless = 0.0;

        double maxFugacityResidual = 0.0;
        double maxPressureResidual = 0.0;
    };

    MaterialList materials_;
    const PengRobinsonEOS& eos_;

    Composition normalizeComposition(
        const Composition& composition
    ) const;

    bool isActive(
        double zi
    ) const;

    Composition wilsonKValues(
        double pressurePa,
        double temperatureK
    ) const;

    double rachfordRiceValue(
        double beta,
        const Composition& z,
        const Composition& k
    ) const;

    double initialBetaEstimate(
        const Composition& z,
        const Composition& k
    ) const;

    GibbsState initialStateFromBeta(
        double pressurePa,
        double temperatureK,
        const Composition& z,
        const Composition& k,
        double betaVapor
    ) const;

    GibbsState initialState(
        double pressurePa,
        double temperatureK,
        const Composition& z
    ) const;

    GibbsEvaluation evaluate(
        const GibbsState& state,
        double pressurePa,
        double temperatureK,
        const Composition& z
    ) const;

    std::vector<double> rhs(
        const GibbsState& state,
        double pressurePa,
        double temperatureK,
        const Composition& z
    ) const;

    RadauStepResult radauIIA3Step(
        const GibbsState& state,
        double pseudoTimeStep,
        double pressurePa,
        double temperatureK,
        const Composition& z
    ) const;

    FlashResult fullRachfordRicePrecheck(
        double pressurePa,
        double temperatureK,
        const Composition& z,
        double tolerance,
        int maxIterations
    ) const;

    double singlePhaseGibbsDimensionless(
        double pressurePa,
        double temperatureK,
        const Composition& z,
        bool vaporPhase
    ) const;

    FlashStatus bestSinglePhase(
        double pressurePa,
        double temperatureK,
        const Composition& z,
        double& bestGibbs
    ) const;

    bool isConverged(
        const GibbsEvaluation& evaluation,
        double tolerance
    ) const;

    FlashResult makeTwoPhaseResult(
        const GibbsEvaluation& evaluation,
        int iterations
    ) const;

    FlashResult makeSinglePhaseResult(
        FlashStatus status,
        double pressurePa,
        double temperatureK,
        const Composition& z
    ) const;

    void fillDensities(
        FlashResult& result,
        double pressurePa,
        double temperatureK
    ) const;

    double averageMolarMassKgPerKmol(
        const Composition& composition
    ) const;
};
