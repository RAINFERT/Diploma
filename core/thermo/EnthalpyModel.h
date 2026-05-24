#pragma once

#include "../Component.h"
#include "Material.h"

#include <array>
#include <cstddef>
#include <string>

struct ComponentEnthalpyData
{
    std::string name;

    double referenceTemperatureK = 298.15;
    double referenceMolarEnthalpyJPerKmol = 0.0;

    // Cp(T) = A + B*T + C*T^2 + D*T^3 + E*T^4
    // J/(kmol*K)
    double cpA = 0.0;
    double cpB = 0.0;
    double cpC = 0.0;
    double cpD = 0.0;
    double cpE = 0.0;
};

using ComponentEnthalpyDataList =
    std::array<ComponentEnthalpyData, ComponentCount>;

ComponentEnthalpyDataList createDefaultEnthalpyData();

class EnthalpyModel
{
public:
    explicit EnthalpyModel(
        const MaterialList& materials
        );

    EnthalpyModel(
        const MaterialList& materials,
        const ComponentEnthalpyDataList& enthalpyData
        );

    const ComponentEnthalpyDataList& data() const;

    ComponentEnthalpyData componentData(
        std::size_t componentIndex
        ) const;

    double componentMolarEnthalpyJPerKmol(
        std::size_t componentIndex,
        double temperatureK
        ) const;

    Composition componentMolarEnthalpiesJPerKmol(
        double temperatureK
        ) const;

    double mixtureMolarEnthalpyJPerKmol(
        double temperatureK,
        const Composition& composition
        ) const;

    double enthalpyFlowW(
        double flowKmolPerS,
        double temperatureK,
        const Composition& composition
        ) const;

    double idealGasMixtureMolarEnthalpyJPerKmol(
        double temperatureK,
        const Composition& composition
        ) const;

    double pengRobinsonDepartureEnthalpyJPerKmol(
        double pressurePa,
        double temperatureK,
        const Composition& composition,
        double compressibilityFactor
        ) const;

    double phaseMolarEnthalpyJPerKmol(
        double pressurePa,
        double temperatureK,
        const Composition& composition,
        double compressibilityFactor
        ) const;

private:
    MaterialList materials_;
    ComponentEnthalpyDataList enthalpyData_;

    void validateTemperature(
        double temperatureK
        ) const;
};
