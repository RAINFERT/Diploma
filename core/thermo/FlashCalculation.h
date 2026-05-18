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

struct FlashResult
{
    FlashStatus status = FlashStatus::NotConverged;

    int iterations = 0;

    // Мольная доля паровой фазы
    double beta = 0.0;

    // Составы фаз
    Composition xLiquid{};
    Composition yVapor{};

    // K_i = y_i / x_i
    Composition kValues{};

    // Z-факторы фаз
    double zLiquid = 0.0;
    double zVapor = 0.0;

    // Коэффициенты летучести
    Composition phiLiquid{};
    Composition phiVapor{};

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
