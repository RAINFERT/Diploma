#pragma once

#include "Material.h"
#include "PengRobinsonEOS.h"

#include <string>

enum class FlashStatus
{
    Converged,
    SingleLiquid,
    SingleVapor,
    NotConverged
};

std::string flashStatusToString(FlashStatus status);

enum class FlashMethod
{
    Unknown,
    RachfordRice,
    GibbsRadau,
    HybridSinglePhasePrecheck,
    HybridFallbackRachfordRice,
    HybridBoundaryRachfordRice,
    HybridFullRachfordRicePrecheck
};

inline const char* flashMethodToString(FlashMethod method)
{
    switch (method)
    {
    case FlashMethod::RachfordRice:
        return "Rachford-Rice K-values";

    case FlashMethod::GibbsRadau:
        return "Gibbs-Radau Gibbs minimization";

    case FlashMethod::HybridSinglePhasePrecheck:
        return "Hybrid: fast single-phase precheck";

    case FlashMethod::HybridFallbackRachfordRice:
        return "Hybrid: fallback Rachford-Rice";

    case FlashMethod::HybridBoundaryRachfordRice:
        return "Hybrid: boundary Rachford-Rice";

    case FlashMethod::HybridFullRachfordRicePrecheck:
        return "Hybrid: full Rachford-Rice precheck";

    case FlashMethod::Unknown:
    default:
        return "Unknown";
    }
}

struct FlashDiagnostics
{
    int pseudoSteps = 0;

    int acceptedRadauSteps = 0;
    int rejectedRadauSteps = 0;

    int totalRadauNewtonIterations = 0;
    int lastRadauNewtonIterations = 0;

    int rhsEvaluations = 0;

    double lastRadauResidualNorm = 0.0;
    double lastRadauStepNorm = 0.0;

    double gibbsTwoPhase = 0.0;
    double gibbsSingleLiquid = 0.0;
    double gibbsSingleVapor = 0.0;

    double maxFugacityResidual = 0.0;
    double maxPressureResidual = 0.0;

    double finalPseudoTimeStep = 0.0;

    double betaLiquidCandidate = 0.0;
    double betaVaporCandidate = 0.0;

    double initialBetaVapor = 0.0;
    double initialGibbsTwoPhase = 0.0;
    double initialMaxFugacityResidual = 0.0;
    double initialMaxPressureResidual = 0.0;

    double precheckBetaVapor = -1.0;
    double precheckF0 = 0.0;
    double precheckF1 = 0.0;

    int precheckStatus = 0;
    int precheckIterations = 0;

    std::string lastRadauFailureMessage;
};


struct FlashResult
{
    FlashStatus status = FlashStatus::NotConverged;

    FlashMethod method = FlashMethod::Unknown;

    int iterations = 0;

    FlashDiagnostics diagnostics{};

    // Мольная доля паровой фазы
    double beta = 0.0;

    // Составы фаз
    Composition xLiquid = makeComposition();
    Composition yVapor = makeComposition();

    Composition kValues = makeComposition();

    Composition phiLiquid = makeComposition();
    Composition phiVapor = makeComposition();

    // Z-факторы фаз
    double zLiquid = 0.0;
    double zVapor = 0.0;

    // Средние молярные массы фаз, kg/kmol
    double molarMassLiquidKgPerKmol = 0.0;
    double molarMassVaporKgPerKmol = 0.0;

    // Молярные плотности фаз, kmol/m3
    double molarDensityLiquidKmolPerM3 = 0.0;
    double molarDensityVaporKmolPerM3 = 0.0;
    double molarDensityMixtureKmolPerM3 = 0.0;

    // Массовые плотности фаз, kg/m3
    double massDensityLiquidKgPerM3 = 0.0;
    double massDensityVaporKgPerM3 = 0.0;
    double massDensityMixtureKgPerM3 = 0.0;

    double molarEnthalpyLiquidJPerKmol = 0.0;
    double molarEnthalpyVaporJPerKmol = 0.0;
    double molarEnthalpyMixtureJPerKmol = 0.0;
};

class FlashCalculation
{
public:
    FlashCalculation(
        const MaterialList& materials,
        const PengRobinsonEOS& eos
        );

    FlashResult calculate(
        double pressurePa,
        double temperatureK,
        const Composition& zOverall,
        double tolerance = 1.0e-8,
        int maxIterations = 100
        ) const;

private:
    static constexpr double R = 8.314462618;

    MaterialList materials_;
    const PengRobinsonEOS& eos_;

    Composition initialKValues(
        double pressurePa,
        double temperatureK,
        const Composition& zOverall
        ) const;

    Composition wilsonKValues(
        double pressurePa,
        double temperatureK
        ) const;

    double rachfordRice(
        double beta,
        const Composition& z,
        const Composition& k
        ) const;

    double rachfordRiceDerivative(
        double beta,
        const Composition& z,
        const Composition& k
        ) const;

    double solveBetaBisection(
        const Composition& z,
        const Composition& k
        ) const;

    double solveBeta(
        const Composition& z,
        const Composition& k
        ) const;

    void computePhaseCompositions(
        double beta,
        const Composition& z,
        const Composition& k,
        Composition& x,
        Composition& y
        ) const;

    double maxRelativeDifference(
        const Composition& a,
        const Composition& b
        ) const;

    double averageMolarMassKgPerKmol(
        const Composition& composition
        ) const;

    double molarDensityKmolPerM3(
        double pressurePa,
        double temperatureK,
        double zFactor
        ) const;

    double massDensityKgPerM3(
        double pressurePa,
        double temperatureK,
        const Composition& composition,
        double zFactor
        ) const;

    bool isAqueousLiquid(
        const Composition& liquidComposition
        ) const;

    double waterMassDensityKgPerM3(
        double pressurePa,
        double temperatureK
        ) const;

    void applyAqueousDensityCorrection(
        FlashResult& result,
        double pressurePa,
        double temperatureK
        ) const;

    void fillDensities(
        FlashResult& result,
        double pressurePa,
        double temperatureK
        ) const;
};
