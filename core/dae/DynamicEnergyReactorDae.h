#pragma once

#include "DynamicStateLayout.h"
#include "../Component.h"
#include "../thermo/DynamicGibbsOdeSystem.h"
#include "../thermo/DynamicEnthalpyModel.h"

#include <vector>

class DynamicWellStirredReactorModel;

using DynamicVector = std::vector<double>;

struct DynamicEnergyReactorDaeEvaluation {
    DynamicComposition massesKg;
    DynamicComposition molesKmol;
    DynamicComposition overallComposition;
    DynamicComposition liquidSplit;
    DynamicGibbsOdeState gibbsOde;

    double totalMolesKmol = 0.0;
    double pressurePa = 0.0;
    double temperatureK = 0.0;
    double energyJ = 0.0;

    double calculatedInventoryEnergyJ = 0.0;

    double calculatedVolumeM3 = 0.0;
    double volumeResidual = 0.0;

    double inventoryEnergyResidual = 0.0;

    double liquidIdealEnthalpyJPerKmol = 0.0;
    double vaporIdealEnthalpyJPerKmol = 0.0;

    double liquidDepartureEnthalpyJPerKmol = 0.0;
    double vaporDepartureEnthalpyJPerKmol = 0.0;

    double liquidTotalEnthalpyJPerKmol = 0.0;
    double vaporTotalEnthalpyJPerKmol = 0.0;

    double twoPhaseMolarEnthalpyJPerKmol = 0.0;
};

class DynamicEnergyReactorDae {
public:
    DynamicEnergyReactorDae(
        const ComponentSet& components,
        const DynamicGibbsOdeSystem& gibbsOde,
        const DynamicEnthalpyModel& enthalpyModel,
        double reactorVolumeM3
        );

    const DynamicStateLayout& layout() const;

    DynamicVector variablesFromState(
        const DynamicComposition& massesKg,
        double energyJ,
        double pressurePa,
        double temperatureK,
        const DynamicComposition& liquidSplit,
        double liquidMolarVolumeM3PerKmol,
        double vaporMolarVolumeM3PerKmol
        ) const;

    DynamicVector zeroVector() const;

    DynamicEnergyReactorDaeEvaluation evaluate(
        const DynamicVector& variables
        ) const;

    DynamicVector computeResidualVector(
        const DynamicVector& variables,
        const DynamicVector& derivatives,
        const DynamicComposition& massRatesKgPerS,
        double energyRateW
        ) const;

    DynamicVector computeCoupledResidualVector(
        const DynamicVector& variables,
        const DynamicVector& derivatives,
        const DynamicWellStirredReactorModel& reactorModel
        ) const;

private:
    const ComponentSet& components_;
    const DynamicGibbsOdeSystem& gibbsOde_;
    const DynamicEnthalpyModel& enthalpyModel_;

    DynamicStateLayout layout_;

    double reactorVolumeM3_;

    void checkVectorSize(
        const DynamicVector& vector,
        const char* vectorName
        ) const;

    void checkCompositionSize(
        const DynamicComposition& composition,
        const char* compositionName
        ) const;

    DynamicComposition extractMassesKg(
        const DynamicVector& variables
        ) const;

    DynamicComposition extractLiquidSplit(
        const DynamicVector& variables
        ) const;

    DynamicComposition massesToMolesKmol(
        const DynamicComposition& massesKg
        ) const;

    DynamicComposition molesToOverallComposition(
        const DynamicComposition& molesKmol
        ) const;

    double extractLiquidMolarVolume(
        const DynamicVector& variables
        ) const;

    double extractVaporMolarVolume(
        const DynamicVector& variables
        ) const;
};
