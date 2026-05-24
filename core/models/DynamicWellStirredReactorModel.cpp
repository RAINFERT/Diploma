#include "DynamicWellStirredReactorModel.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>

DynamicWellStirredReactorModel::DynamicWellStirredReactorModel(
    const ComponentSet& components,
    const DynamicEnthalpyModel& enthalpyModel,
    DynamicWellStirredReactorParameters parameters
    )
    : components_(components),
    enthalpyModel_(enthalpyModel),
    parameters_(std::move(parameters))
{
    if (components_.size() == 0) {
        throw std::runtime_error(
            "DynamicWellStirredReactorModel: empty component set"
            );
    }

    if (parameters_.volumeM3 <= 0.0) {
        throw std::runtime_error(
            "DynamicWellStirredReactorModel: reactor volume must be positive"
            );
    }

    validateComposition(
        parameters_.inletComposition,
        "inletComposition"
        );
}

const DynamicWellStirredReactorParameters&
DynamicWellStirredReactorModel::parameters() const {
    return parameters_;
}

void DynamicWellStirredReactorModel::validateComposition(
    const DynamicComposition& composition,
    const char* name
    ) const {
    if (composition.size() != components_.size()) {
        throw std::runtime_error(
            std::string("DynamicWellStirredReactorModel: invalid size of ") +
            name
            );
    }

    double sum = 0.0;

    for (double value : composition) {
        if (value < 0.0) {
            throw std::runtime_error(
                std::string("DynamicWellStirredReactorModel: negative value in ") +
                name
                );
        }

        sum += value;
    }

    if (sum <= 0.0) {
        throw std::runtime_error(
            std::string("DynamicWellStirredReactorModel: zero sum in ") +
            name
            );
    }
}

DynamicMassBalanceResult
DynamicWellStirredReactorModel::computeMassBalance(
    double pressurePa,
    const DynamicComposition& massesKg,
    const DynamicComposition& vaporComposition
    ) const {
    if (massesKg.size() != components_.size()) {
        throw std::runtime_error(
            "DynamicWellStirredReactorModel: invalid massesKg size"
            );
    }

    validateComposition(
        vaporComposition,
        "vaporComposition"
        );

    DynamicMassBalanceResult result;
    result.massRatesKgPerS.assign(
        components_.size(),
        0.0
        );

    double outletFlowKmolPerS =
        parameters_.outletValveCoefficientKmolPerSPa
        * (pressurePa - parameters_.inletPressurePa);

    if (!parameters_.allowReverseOutletFlow) {
        outletFlowKmolPerS =
            std::max(0.0, outletFlowKmolPerS);
    }

    result.outletMolarFlowKmolPerS =
        outletFlowKmolPerS;

    for (std::size_t i = 0; i < components_.size(); ++i) {
        const double molecularWeight =
            components_[i].molarMassKgPerKmol;

        const double inletMassRate =
            parameters_.inletMolarFlowKmolPerS
            * parameters_.inletComposition[i]
            * molecularWeight;

        const double outletMassRate =
            outletFlowKmolPerS
            * vaporComposition[i]
            * molecularWeight;

        result.massRatesKgPerS[i] =
            inletMassRate - outletMassRate;
    }

    return result;
}

DynamicEnergyBalanceResult
DynamicWellStirredReactorModel::computeEnergyBalance(
    double reactorPressurePa,
    double reactorTemperatureK,
    double outletMolarFlowKmolPerS,
    const DynamicComposition& outletComposition
    ) const {
    if (reactorPressurePa <= 0.0) {
        throw std::runtime_error(
            "DynamicWellStirredReactorModel: reactor pressure must be positive"
            );
    }

    if (reactorTemperatureK <= 0.0) {
        throw std::runtime_error(
            "DynamicWellStirredReactorModel: reactor temperature must be positive"
            );
    }

    validateComposition(
        outletComposition,
        "outletComposition"
        );

    DynamicEnergyBalanceResult result;

    // Вход пока считаем паровой фазой.
    // Теперь здесь учитывается PR departure через DynamicEnthalpyModel.
    result.inletEnthalpyFlowW =
        enthalpyModel_.enthalpyFlowW(
            parameters_.inletMolarFlowKmolPerS,
            parameters_.inletPressurePa,
            parameters_.inletTemperatureK,
            parameters_.inletComposition,
            DynamicPhase::Vapor
            );

    // Выход берем из паровой фазы реактора.
    // Теперь Hdot_out тоже учитывает departure.
    result.outletEnthalpyFlowW =
        enthalpyModel_.enthalpyFlowW(
            outletMolarFlowKmolPerS,
            reactorPressurePa,
            reactorTemperatureK,
            outletComposition,
            DynamicPhase::Vapor
            );

    result.heatTransferW =
        parameters_.heatTransferCoefficientWPerM2K
        * parameters_.heatTransferAreaM2
        * (parameters_.jacketTemperatureK - reactorTemperatureK);

    result.energyRateW =
        result.heatTransferW
        + result.inletEnthalpyFlowW
        - result.outletEnthalpyFlowW;

    return result;
}
