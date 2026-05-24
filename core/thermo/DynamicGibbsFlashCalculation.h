#pragma once

#include "DynamicPengRobinsonEOS.h"

#include <string>
#include <vector>

struct DynamicGibbsFlashOptions {
    double initialLiquidSplit = 0.5;

    double pseudoTimeStep = 0.05;
    double mobility = 1.0;

    double tolerance = 1.0e-8;
    int maxIterations = 20000;

    double minSplit = 1.0e-12;
    double minPhaseFraction = 1.0e-10;
    double maxSplitStep = 0.05;

    bool useGibbsProbabilityFactor = false;
};

struct DynamicGibbsFlashResult {
    bool converged = false;
    bool singlePhase = false;

    std::string status;

    int iterations = 0;

    double betaLiquid = 0.0;
    double betaVapor = 0.0;

    DynamicComposition liquidSplit;

    DynamicComposition xLiquid;
    DynamicComposition yVapor;

    DynamicComposition phiLiquid;
    DynamicComposition phiVapor;

    DynamicZFactors zLiquid;
    DynamicZFactors zVapor;

    double maxFugacityResidual = 0.0;
    double gibbsDeviationDimensionless = 0.0;
};

class DynamicGibbsFlashCalculation {
public:
    explicit DynamicGibbsFlashCalculation(
        const DynamicPengRobinsonEOS& eos
        );

    DynamicGibbsFlashResult calculate(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& overallComposition,
        const DynamicGibbsFlashOptions& options = DynamicGibbsFlashOptions{}
        ) const;

    DynamicGibbsFlashResult evaluateEquilibriumResidual(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& overallComposition,
        const DynamicComposition& liquidSplit
        ) const;

    DynamicComposition fugacityResiduals(
        double pressurePa,
        const DynamicComposition& xLiquid,
        const DynamicComposition& yVapor,
        const DynamicComposition& phiLiquid,
        const DynamicComposition& phiVapor
        ) const;

private:
    const DynamicPengRobinsonEOS& eos_;

    DynamicComposition normalizeComposition(
        const DynamicComposition& composition
        ) const;

    DynamicGibbsFlashResult evaluateState(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& overallComposition,
        const DynamicComposition& liquidSplit
        ) const;

    DynamicComposition computeFugacityResiduals(
        double pressurePa,
        const DynamicComposition& xLiquid,
        const DynamicComposition& yVapor,
        const DynamicComposition& phiLiquid,
        const DynamicComposition& phiVapor
        ) const;

    double computeGibbsDeviationDimensionless(
        const DynamicComposition& overallComposition,
        const DynamicComposition& liquidSplit,
        const DynamicGibbsFlashResult& state
        ) const;
};
