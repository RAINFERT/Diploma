#pragma once

#include "../Component.h"
#include "DynamicPengRobinsonEOS.h"

#include <cstddef>
#include <string>

enum class DynamicPhase {
    Vapor,
    Liquid
};

struct DynamicComponentEnthalpyInfo {
    double vaporCpJPerKmolK = 0.0;
    double liquidCpJPerKmolK = 0.0;

    double vaporH_JPerKmol = 0.0;
    double liquidH_JPerKmol = 0.0;

    std::string vaporSource;
    std::string liquidSource;
};

class DynamicEnthalpyModel {
public:
    explicit DynamicEnthalpyModel(
        const ComponentSet& components
        );

    DynamicEnthalpyModel(
        const ComponentSet& components,
        const DynamicPengRobinsonEOS& eos
        );

    double referenceTemperatureK() const;

    DynamicComponentEnthalpyInfo componentInfo(
        std::size_t componentIndex,
        double temperatureK
        ) const;

    double componentVaporHeatCapacityJPerKmolK(
        std::size_t componentIndex,
        double temperatureK
        ) const;

    double componentLiquidHeatCapacityJPerKmolK(
        std::size_t componentIndex,
        double temperatureK
        ) const;

    double componentVaporEnthalpyJPerKmol(
        std::size_t componentIndex,
        double temperatureK
        ) const;

    double componentLiquidEnthalpyJPerKmol(
        std::size_t componentIndex,
        double temperatureK
        ) const;

    double vaporMixtureMolarEnthalpyJPerKmol(
        double temperatureK,
        const DynamicComposition& yVapor
        ) const;

    double liquidMixtureMolarEnthalpyJPerKmol(
        double temperatureK,
        const DynamicComposition& xLiquid
        ) const;

    double mixtureMolarEnthalpyJPerKmol(
        double temperatureK,
        const DynamicComposition& composition
        ) const;

    double twoPhaseMolarEnthalpyJPerKmol(
        double temperatureK,
        double betaLiquid,
        double betaVapor,
        const DynamicComposition& xLiquid,
        const DynamicComposition& yVapor
        ) const;

    double inventoryEnergyJ(
        double totalMolesKmol,
        double temperatureK,
        double betaLiquid,
        double betaVapor,
        const DynamicComposition& xLiquid,
        const DynamicComposition& yVapor
        ) const;

    double enthalpyFlowW(
        double molarFlowKmolPerS,
        double temperatureK,
        const DynamicComposition& composition
        ) const;

    double phaseMolarEnthalpyJPerKmol(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& composition,
        DynamicPhase phase
        ) const;

    double enthalpyFlowW(
        double molarFlowKmolPerS,
        double pressurePa,
        double temperatureK,
        const DynamicComposition& composition,
        DynamicPhase phase
        ) const;

    double phaseMolarEnthalpyJPerKmol(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& composition,
        double zFactor
        ) const;

    double phaseDepartureEnthalpyJPerKmol(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& composition,
        double zFactor
        ) const;

    double twoPhaseMolarEnthalpyJPerKmol(
        double pressurePa,
        double temperatureK,
        double betaLiquid,
        double betaVapor,
        const DynamicComposition& xLiquid,
        const DynamicComposition& yVapor,
        double zLiquid,
        double zVapor
        ) const;

    double inventoryEnergyJ(
        double totalMolesKmol,
        double pressurePa,
        double temperatureK,
        double betaLiquid,
        double betaVapor,
        const DynamicComposition& xLiquid,
        const DynamicComposition& yVapor,
        double zLiquid,
        double zVapor
        ) const;

private:
    const ComponentSet& components_;

    static constexpr double ReferenceTemperatureK = 298.15;

    const DynamicPengRobinsonEOS* eos_ = nullptr;

    void checkComponentIndex(
        std::size_t componentIndex
        ) const;

    void checkComposition(
        const DynamicComposition& composition
        ) const;

    bool hasVaporCpCorrelation(
        const ComponentProperties& component
        ) const;

    bool hasLiquidCpCorrelation(
        const ComponentProperties& component
        ) const;

    double evaluatePolynomial100(
        double A,
        double B,
        double C,
        double D,
        double E,
        double temperatureK
        ) const;

    double fallbackVaporHeatCapacityJPerKmolK(
        const ComponentProperties& component
        ) const;

    double fallbackLiquidHeatCapacityJPerKmolK(
        const ComponentProperties& component
        ) const;

    double integrateComponentCpJPerKmol(
        std::size_t componentIndex,
        double temperatureK,
        DynamicPhase phase
        ) const;
};
