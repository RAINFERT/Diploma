#include "DynamicEnthalpyModel.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

DynamicEnthalpyModel::DynamicEnthalpyModel(
    const ComponentSet& components
    )
    : components_(components)
{
    if (components_.size() == 0) {
        throw std::runtime_error(
            "DynamicEnthalpyModel: empty component set"
            );
    }
}

DynamicEnthalpyModel::DynamicEnthalpyModel(
    const ComponentSet& components,
    const DynamicPengRobinsonEOS& eos
    )
    : components_(components),
    eos_(&eos)
{
    if (components_.size() == 0) {
        throw std::runtime_error(
            "DynamicEnthalpyModel: empty component set"
            );
    }
}

double DynamicEnthalpyModel::referenceTemperatureK() const {
    return ReferenceTemperatureK;
}

void DynamicEnthalpyModel::checkComponentIndex(
    std::size_t componentIndex
    ) const {
    if (componentIndex >= components_.size()) {
        throw std::runtime_error(
            "DynamicEnthalpyModel: component index out of range"
            );
    }
}

void DynamicEnthalpyModel::checkComposition(
    const DynamicComposition& composition
    ) const {
    if (composition.size() != components_.size()) {
        throw std::runtime_error(
            "DynamicEnthalpyModel: composition size mismatch"
            );
    }

    double sum = 0.0;

    for (double value : composition) {
        if (value < 0.0) {
            throw std::runtime_error(
                "DynamicEnthalpyModel: negative composition value"
                );
        }

        sum += value;
    }

    if (sum <= 0.0) {
        throw std::runtime_error(
            "DynamicEnthalpyModel: composition sum must be positive"
            );
    }
}

double DynamicEnthalpyModel::evaluatePolynomial100(
    double A,
    double B,
    double C,
    double D,
    double E,
    double temperatureK
    ) const {
    const double T = temperatureK;

    return A
           + B * T
           + C * T * T
           + D * T * T * T
           + E * T * T * T * T;
}

bool DynamicEnthalpyModel::hasVaporCpCorrelation(
    const ComponentProperties& component
    ) const {
    const double magnitude =
        std::abs(component.idealGasCpA)
        + std::abs(component.idealGasCpB)
        + std::abs(component.idealGasCpC)
        + std::abs(component.idealGasCpD)
        + std::abs(component.idealGasCpE);

    return component.idealGasCpEqno == 100
           && magnitude > 0.0;
}

bool DynamicEnthalpyModel::hasLiquidCpCorrelation(
    const ComponentProperties& component
    ) const {
    const double magnitude =
        std::abs(component.liquidCpA)
        + std::abs(component.liquidCpB)
        + std::abs(component.liquidCpC)
        + std::abs(component.liquidCpD)
        + std::abs(component.liquidCpE);

    // Пока безопасно используем только eqno = 100.
    // Для eqno = 16 нужна отдельная формула; угадывать ее нельзя.
    return component.liquidCpEqno == 100
           && magnitude > 0.0;
}

double DynamicEnthalpyModel::fallbackVaporHeatCapacityJPerKmolK(
    const ComponentProperties& component
    ) const {
    if (component.name == "Ethane") {
        return 52000.0;
    }

    if (component.name == "N-pentane") {
        return 120000.0;
    }

    if (component.name == "Water") {
        return 33600.0;
    }

    return 100000.0;
}

double DynamicEnthalpyModel::fallbackLiquidHeatCapacityJPerKmolK(
    const ComponentProperties& component
    ) const {
    if (component.name == "Ethane") {
        return 110000.0;
    }

    if (component.name == "N-pentane") {
        return 167000.0;
    }

    if (component.name == "Water") {
        return 75300.0;
    }

    return 120000.0;
}

double DynamicEnthalpyModel::componentVaporHeatCapacityJPerKmolK(
    std::size_t componentIndex,
    double temperatureK
    ) const {
    checkComponentIndex(componentIndex);

    if (temperatureK <= 0.0) {
        throw std::runtime_error(
            "DynamicEnthalpyModel: temperature must be positive"
            );
    }

    const ComponentProperties& component =
        components_[componentIndex];

    double cp = 0.0;

    if (hasVaporCpCorrelation(component)) {
        cp =
            evaluatePolynomial100(
                component.idealGasCpA,
                component.idealGasCpB,
                component.idealGasCpC,
                component.idealGasCpD,
                component.idealGasCpE,
                temperatureK
                );
    } else {
        cp =
            fallbackVaporHeatCapacityJPerKmolK(
                component
                );
    }

    if (!std::isfinite(cp) || cp <= 1000.0 || cp > 1.0e7) {
        return fallbackVaporHeatCapacityJPerKmolK(
            component
            );
    }

    return cp;
}

double DynamicEnthalpyModel::componentLiquidHeatCapacityJPerKmolK(
    std::size_t componentIndex,
    double temperatureK
    ) const {
    checkComponentIndex(componentIndex);

    if (temperatureK <= 0.0) {
        throw std::runtime_error(
            "DynamicEnthalpyModel: temperature must be positive"
            );
    }

    const ComponentProperties& component =
        components_[componentIndex];

    double cp = 0.0;

    if (hasLiquidCpCorrelation(component)) {
        cp =
            evaluatePolynomial100(
                component.liquidCpA,
                component.liquidCpB,
                component.liquidCpC,
                component.liquidCpD,
                component.liquidCpE,
                temperatureK
                );
    } else {
        cp =
            fallbackLiquidHeatCapacityJPerKmolK(
                component
                );
    }

    if (!std::isfinite(cp) || cp <= 1000.0 || cp > 1.0e7) {
        return fallbackLiquidHeatCapacityJPerKmolK(
            component
            );
    }

    return cp;
}

double DynamicEnthalpyModel::integrateComponentCpJPerKmol(
    std::size_t componentIndex,
    double temperatureK,
    DynamicPhase phase
    ) const {
    checkComponentIndex(componentIndex);

    if (temperatureK <= 0.0) {
        throw std::runtime_error(
            "DynamicEnthalpyModel: temperature must be positive"
            );
    }

    const double T0 = ReferenceTemperatureK;
    const double T1 = temperatureK;

    if (std::abs(T1 - T0) < 1.0e-12) {
        return 0.0;
    }

    constexpr int SegmentCount = 64;

    const double a = std::min(T0, T1);
    const double b = std::max(T0, T1);
    const double h = (b - a) / SegmentCount;

    double integral = 0.0;

    for (int k = 0; k <= SegmentCount; ++k) {
        const double T = a + k * h;

        double cp = 0.0;

        if (phase == DynamicPhase::Vapor) {
            cp =
                componentVaporHeatCapacityJPerKmolK(
                    componentIndex,
                    T
                    );
        } else {
            cp =
                componentLiquidHeatCapacityJPerKmolK(
                    componentIndex,
                    T
                    );
        }

        double weight = 2.0;

        if (k == 0 || k == SegmentCount) {
            weight = 1.0;
        } else if (k % 2 == 1) {
            weight = 4.0;
        }

        integral += weight * cp;
    }

    integral *= h / 3.0;

    if (T1 < T0) {
        integral = -integral;
    }

    return integral;
}

double DynamicEnthalpyModel::componentVaporEnthalpyJPerKmol(
    std::size_t componentIndex,
    double temperatureK
    ) const {
    return integrateComponentCpJPerKmol(
        componentIndex,
        temperatureK,
        DynamicPhase::Vapor
        );
}

double DynamicEnthalpyModel::componentLiquidEnthalpyJPerKmol(
    std::size_t componentIndex,
    double temperatureK
    ) const {
    return integrateComponentCpJPerKmol(
        componentIndex,
        temperatureK,
        DynamicPhase::Liquid
        );
}

DynamicComponentEnthalpyInfo DynamicEnthalpyModel::componentInfo(
    std::size_t componentIndex,
    double temperatureK
    ) const {
    checkComponentIndex(componentIndex);

    const ComponentProperties& component =
        components_[componentIndex];

    DynamicComponentEnthalpyInfo info;

    info.vaporCpJPerKmolK =
        componentVaporHeatCapacityJPerKmolK(
            componentIndex,
            temperatureK
            );

    info.liquidCpJPerKmolK =
        componentLiquidHeatCapacityJPerKmolK(
            componentIndex,
            temperatureK
            );

    info.vaporH_JPerKmol =
        componentVaporEnthalpyJPerKmol(
            componentIndex,
            temperatureK
            );

    info.liquidH_JPerKmol =
        componentLiquidEnthalpyJPerKmol(
            componentIndex,
            temperatureK
            );

    info.vaporSource =
        hasVaporCpCorrelation(component)
            ? "loaded vapor ideal-gas Cp polynomial"
            : "fallback vapor Cp";

    info.liquidSource =
        hasLiquidCpCorrelation(component)
            ? "loaded liquid Cp polynomial"
            : "fallback liquid Cp";

    return info;
}

double DynamicEnthalpyModel::vaporMixtureMolarEnthalpyJPerKmol(
    double temperatureK,
    const DynamicComposition& yVapor
    ) const {
    checkComposition(yVapor);

    const double sum =
        std::accumulate(
            yVapor.begin(),
            yVapor.end(),
            0.0
            );

    double h = 0.0;

    for (std::size_t i = 0; i < components_.size(); ++i) {
        h +=
            yVapor[i] / sum
            * componentVaporEnthalpyJPerKmol(
                i,
                temperatureK
                );
    }

    return h;
}

double DynamicEnthalpyModel::liquidMixtureMolarEnthalpyJPerKmol(
    double temperatureK,
    const DynamicComposition& xLiquid
    ) const {
    checkComposition(xLiquid);

    const double sum =
        std::accumulate(
            xLiquid.begin(),
            xLiquid.end(),
            0.0
            );

    double h = 0.0;

    for (std::size_t i = 0; i < components_.size(); ++i) {
        h +=
            xLiquid[i] / sum
            * componentLiquidEnthalpyJPerKmol(
                i,
                temperatureK
                );
    }

    return h;
}

double DynamicEnthalpyModel::mixtureMolarEnthalpyJPerKmol(
    double temperatureK,
    const DynamicComposition& composition
    ) const {
    // Для однофазного входного потока пока считаем его vapor/ideal-gas
    // энтальпией. Позже можно добавить признак phase для feed.
    return vaporMixtureMolarEnthalpyJPerKmol(
        temperatureK,
        composition
        );
}

double DynamicEnthalpyModel::twoPhaseMolarEnthalpyJPerKmol(
    double temperatureK,
    double betaLiquid,
    double betaVapor,
    const DynamicComposition& xLiquid,
    const DynamicComposition& yVapor
    ) const {
    const double hLiquid =
        liquidMixtureMolarEnthalpyJPerKmol(
            temperatureK,
            xLiquid
            );

    const double hVapor =
        vaporMixtureMolarEnthalpyJPerKmol(
            temperatureK,
            yVapor
            );

    return betaLiquid * hLiquid
           + betaVapor * hVapor;
}

double DynamicEnthalpyModel::inventoryEnergyJ(
    double totalMolesKmol,
    double temperatureK,
    double betaLiquid,
    double betaVapor,
    const DynamicComposition& xLiquid,
    const DynamicComposition& yVapor
    ) const {
    const double hMix =
        twoPhaseMolarEnthalpyJPerKmol(
            temperatureK,
            betaLiquid,
            betaVapor,
            xLiquid,
            yVapor
            );

    return totalMolesKmol * hMix;
}

double DynamicEnthalpyModel::enthalpyFlowW(
    double molarFlowKmolPerS,
    double temperatureK,
    const DynamicComposition& composition
    ) const {
    return molarFlowKmolPerS
           * mixtureMolarEnthalpyJPerKmol(
               temperatureK,
               composition
               );
}

double DynamicEnthalpyModel::phaseDepartureEnthalpyJPerKmol(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& composition,
    double zFactor
    ) const {
    if (eos_ == nullptr) {
        return 0.0;
    }

    return eos_->computeDepartureEnthalpyJPerKmol(
        pressurePa,
        temperatureK,
        composition,
        zFactor
        );
}

double DynamicEnthalpyModel::phaseMolarEnthalpyJPerKmol(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& composition,
    double zFactor
    ) const {
    const double idealPart =
        vaporMixtureMolarEnthalpyJPerKmol(
            temperatureK,
            composition
            );

    const double departurePart =
        phaseDepartureEnthalpyJPerKmol(
            pressurePa,
            temperatureK,
            composition,
            zFactor
            );

    return idealPart + departurePart;
}

double DynamicEnthalpyModel::twoPhaseMolarEnthalpyJPerKmol(
    double pressurePa,
    double temperatureK,
    double betaLiquid,
    double betaVapor,
    const DynamicComposition& xLiquid,
    const DynamicComposition& yVapor,
    double zLiquid,
    double zVapor
    ) const {
    const double hLiquid =
        phaseMolarEnthalpyJPerKmol(
            pressurePa,
            temperatureK,
            xLiquid,
            zLiquid
            );

    const double hVapor =
        phaseMolarEnthalpyJPerKmol(
            pressurePa,
            temperatureK,
            yVapor,
            zVapor
            );

    return betaLiquid * hLiquid
           + betaVapor * hVapor;
}

double DynamicEnthalpyModel::inventoryEnergyJ(
    double totalMolesKmol,
    double pressurePa,
    double temperatureK,
    double betaLiquid,
    double betaVapor,
    const DynamicComposition& xLiquid,
    const DynamicComposition& yVapor,
    double zLiquid,
    double zVapor
    ) const {
    const double hMix =
        twoPhaseMolarEnthalpyJPerKmol(
            pressurePa,
            temperatureK,
            betaLiquid,
            betaVapor,
            xLiquid,
            yVapor,
            zLiquid,
            zVapor
            );

    return totalMolesKmol * hMix;
}

double DynamicEnthalpyModel::phaseMolarEnthalpyJPerKmol(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& composition,
    DynamicPhase phase
    ) const {
    if (eos_ == nullptr) {
        if (phase == DynamicPhase::Vapor) {
            return vaporMixtureMolarEnthalpyJPerKmol(
                temperatureK,
                composition
                );
        }

        return liquidMixtureMolarEnthalpyJPerKmol(
            temperatureK,
            composition
            );
    }

    const DynamicZFactors zFactors =
        eos_->computeZFactors(
            pressurePa,
            temperatureK,
            composition
            );

    const double zFactor =
        phase == DynamicPhase::Vapor
            ? zFactors.vapor
            : zFactors.liquid;

    return phaseMolarEnthalpyJPerKmol(
        pressurePa,
        temperatureK,
        composition,
        zFactor
        );
}

double DynamicEnthalpyModel::enthalpyFlowW(
    double molarFlowKmolPerS,
    double pressurePa,
    double temperatureK,
    const DynamicComposition& composition,
    DynamicPhase phase
    ) const {
    return molarFlowKmolPerS
           * phaseMolarEnthalpyJPerKmol(
               pressurePa,
               temperatureK,
               composition,
               phase
               );
}
