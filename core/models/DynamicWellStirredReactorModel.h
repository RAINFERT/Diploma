#pragma once

#include "../Component.h"
#include "../thermo/DynamicGibbsFlashCalculation.h"
#include "../thermo/DynamicEnthalpyModel.h"

#include <vector>

struct DynamicWellStirredReactorParameters {
    double volumeM3 = 1.0;

    double inletMolarFlowKmolPerS = 0.0;
    double outletValveCoefficientKmolPerSPa = 0.0;

    double inletPressurePa = 101325.0;
    double inletTemperatureK = 298.15;

    DynamicComposition inletComposition;

    bool allowReverseOutletFlow = false;

    double heatTransferAreaM2 = 0.0;
    double heatTransferCoefficientWPerM2K = 0.0;
    double jacketTemperatureK = 298.15;

};

struct DynamicMassBalanceResult {
    DynamicComposition massRatesKgPerS;
    double outletMolarFlowKmolPerS = 0.0;
};

struct DynamicEnergyBalanceResult {
    double energyRateW = 0.0;

    double heatTransferW = 0.0;
    double inletEnthalpyFlowW = 0.0;
    double outletEnthalpyFlowW = 0.0;
};

class DynamicWellStirredReactorModel {
public:
    DynamicWellStirredReactorModel(
        const ComponentSet& components,
        const DynamicEnthalpyModel& enthalpyModel,
        DynamicWellStirredReactorParameters parameters
        );

    const DynamicWellStirredReactorParameters& parameters() const;

    DynamicMassBalanceResult computeMassBalance(
        double pressurePa,
        const DynamicComposition& massesKg,
        const DynamicComposition& vaporComposition
        ) const;

    DynamicEnergyBalanceResult computeEnergyBalance(
        double reactorPressurePa,
        double reactorTemperatureK,
        double outletMolarFlowKmolPerS,
        const DynamicComposition& outletComposition
        ) const;

private:
    const ComponentSet& components_;
    DynamicWellStirredReactorParameters parameters_;

    void validateComposition(
        const DynamicComposition& composition,
        const char* name
        ) const;

    const DynamicEnthalpyModel& enthalpyModel_;
};
