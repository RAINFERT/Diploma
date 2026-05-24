#pragma once

#include "DynamicPengRobinsonEOS.h"

#include <vector>

struct DynamicGibbsOdeOptions {
    double minSplit = 1.0e-12;
    double minPhaseFraction = 1.0e-12;

    double compositionMobility = 1.0;
    double volumeMobility = 1.0;

    double maxExpArgument = 40.0;
    double maxPressureResidual = 100.0;
};

struct DynamicGibbsOdeState {
    DynamicComposition overallComposition;

    DynamicComposition liquidSplit;

    DynamicComposition liquidMolesFraction;
    DynamicComposition vaporMolesFraction;

    DynamicComposition xLiquid;
    DynamicComposition yVapor;

    DynamicComposition phiLiquid;
    DynamicComposition phiVapor;

    DynamicComposition fugacityLogResiduals;

    // Right-hand side of article ODE for s_i = N_i^L / N_i^0
    DynamicComposition liquidSplitRates;

    double betaLiquid = 0.0;
    double betaVapor = 0.0;

    double liquidMolarVolumeM3PerKmol = 0.0;
    double vaporMolarVolumeM3PerKmol = 0.0;

    double liquidPressureFromEosPa = 0.0;
    double vaporPressureFromEosPa = 0.0;

    double liquidZ = 0.0;
    double vaporZ = 0.0;

    // Right-hand side of article ODE for molar volumes
    double liquidMolarVolumeRate = 0.0;
    double vaporMolarVolumeRate = 0.0;

    double gibbsDeviationDimensionless = 0.0;
    double maxFugacityResidual = 0.0;
};

class DynamicGibbsOdeSystem {
public:
    explicit DynamicGibbsOdeSystem(
        const DynamicPengRobinsonEOS& eos
        );

    DynamicGibbsOdeState evaluate(
        double pressurePa,
        double temperatureK,
        double totalMolesKmol,
        const DynamicComposition& overallComposition,
        const DynamicComposition& liquidSplit,
        double liquidMolarVolumeM3PerKmol,
        double vaporMolarVolumeM3PerKmol,
        const DynamicGibbsOdeOptions& options = DynamicGibbsOdeOptions{}
        ) const;

private:
    const DynamicPengRobinsonEOS& eos_;

    DynamicComposition normalizeComposition(
        const DynamicComposition& composition
        ) const;

    double safePositive(
        double value
        ) const;

    double clampValue(
        double value,
        double minValue,
        double maxValue
        ) const;

    double pseudoCriticalVolumeM3PerKmol(
        const DynamicComposition& composition
        ) const;
};
