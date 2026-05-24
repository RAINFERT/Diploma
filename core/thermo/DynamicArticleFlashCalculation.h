#pragma once

#include "DynamicGibbsOdeSystem.h"
#include "DynamicPengRobinsonEOS.h"

#include <string>

struct DynamicArticleFlashOptions {
    double initialLiquidSplit = 0.5;

    double pseudoTimeStep = 0.02;
    int maxIterations = 50000;

    double fugacityTolerance = 1.0e-8;
    double pressureTolerance = 1.0e-8;
    double rateTolerance = 1.0e-10;

    double minSplit = 1.0e-12;
    double singlePhaseTolerance = 1.0e-8;

    double compositionMobility = 1.0;
    double volumeMobility = 1.0;
};

enum class DynamicArticlePhaseState {
    TwoPhase,
    SingleLiquid,
    SingleVapor,
    NotConverged
};

struct DynamicArticleFlashResult {
    bool converged = false;
    std::string status;

    DynamicArticlePhaseState phaseState =
        DynamicArticlePhaseState::NotConverged;

    int iterations = 0;

    double betaLiquid = 0.0;
    double betaVapor = 0.0;

    DynamicComposition liquidSplit;

    DynamicComposition xLiquid;
    DynamicComposition yVapor;

    DynamicComposition phiLiquid;
    DynamicComposition phiVapor;

    DynamicComposition fugacityLogResiduals;

    double liquidMolarVolumeM3PerKmol = 0.0;
    double vaporMolarVolumeM3PerKmol = 0.0;

    double liquidZ = 0.0;
    double vaporZ = 0.0;

    double liquidPressureFromEosPa = 0.0;
    double vaporPressureFromEosPa = 0.0;

    double gibbsDeviationDimensionless = 0.0;
    double maxFugacityResidual = 0.0;
    double maxPressureResidual = 0.0;
    double maxRateResidual = 0.0;
};

class DynamicArticleFlashCalculation {
public:
    explicit DynamicArticleFlashCalculation(
        const DynamicPengRobinsonEOS& eos
        );

    DynamicArticleFlashResult calculate(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& overallComposition,
        const DynamicArticleFlashOptions& options =
        DynamicArticleFlashOptions{}
        ) const;

private:
    const DynamicPengRobinsonEOS& eos_;
    DynamicGibbsOdeSystem odeSystem_;

    DynamicComposition normalizeComposition(
        const DynamicComposition& composition
        ) const;

    DynamicComposition initialLiquidSplit(
        const DynamicComposition& z,
        const DynamicArticleFlashOptions& options
        ) const;

    void initialMolarVolumes(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& z,
        double& liquidMolarVolumeM3PerKmol,
        double& vaporMolarVolumeM3PerKmol
        ) const;

    DynamicArticleFlashResult buildResult(
        double pressurePa,
        double temperatureK,
        const DynamicGibbsOdeState& state,
        int iterations,
        bool converged,
        const std::string& status,
        const DynamicArticleFlashOptions& options
        ) const;

    bool isConvergedTwoPhase(
        const DynamicGibbsOdeState& state,
        double pressurePa,
        const DynamicArticleFlashOptions& options,
        double& maxPressureResidual,
        double& maxRateResidual
        ) const;

    bool isSingleLiquidBoundary(
        const DynamicGibbsOdeState& state,
        const DynamicArticleFlashOptions& options
        ) const;

    bool isSingleVaporBoundary(
        const DynamicGibbsOdeState& state,
        const DynamicArticleFlashOptions& options
        ) const;

    double clampValue(
        double value,
        double minValue,
        double maxValue
        ) const;
};
