#include "ChemsepMaterialBuilder.h"

#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace {

double requireFiniteScalar(
    const ChemsepComponent& component,
    const std::string& propertyName
    )
{
    const double value =
        component.scalar(propertyName);

    if (!std::isfinite(value)) {
        throw std::runtime_error(
            "ChemSep scalar property is not finite: " +
            component.componentName() +
            "." +
            propertyName
            );
    }

    return value;
}

const ChemsepCorrelation& requireCorrelation(
    const ChemsepComponent& component,
    const std::string& correlationName
    )
{
    if (!component.hasCorrelation(correlationName)) {
        throw std::runtime_error(
            "ChemSep correlation is missing: " +
            component.componentName() +
            "." +
            correlationName
            );
    }

    return component.correlation(correlationName);
}

AntoineCoefficients makeVapourPressureCoefficients(
    const ChemsepComponent& component
    )
{
    const ChemsepCorrelation& vapourPressure =
        requireCorrelation(
            component,
            "vapour_pressure"
            );

    if (vapourPressure.eqno != 101) {
        std::ostringstream message;
        message
            << "Unsupported ChemSep vapour_pressure_eqno for "
            << component.componentName()
            << ". Expected 101, got "
            << vapourPressure.eqno;

        throw std::runtime_error(message.str());
    }

    if (
        !std::isfinite(vapourPressure.A) ||
        !std::isfinite(vapourPressure.B) ||
        !std::isfinite(vapourPressure.C) ||
        !std::isfinite(vapourPressure.D) ||
        !std::isfinite(vapourPressure.E)
        ) {
        throw std::runtime_error(
            "ChemSep vapour_pressure coefficients are incomplete for " +
            component.componentName()
            );
    }

    // Текущая структура AntoineCoefficients в Material.h использует:
    //
    // ln(P[Pa]) = A + B / (C + T) + D * ln(T) + E * T^F
    //
    // Для ChemSep vapour_pressure_eqno = 101:
    //
    // ln(P[Pa]) = A + B / T + C_chemsep * ln(T) + D_chemsep * T^E_chemsep
    //
    // Поэтому:
    // A = vapour_pressure_A
    // B = vapour_pressure_B
    // C = 0.0
    // D = vapour_pressure_C
    // E = vapour_pressure_D
    // F = vapour_pressure_E

    return AntoineCoefficients{
        vapourPressure.A,
        vapourPressure.B,
        0.0,
        vapourPressure.C,
        vapourPressure.D,
        vapourPressure.E
    };
}

void validateMaterial(
    const Material& material
    )
{
    if (material.criticalTemperatureK <= 0.0) {
        throw std::runtime_error(
            "Invalid critical temperature for material: " +
            material.name
            );
    }

    if (material.criticalPressurePa <= 0.0) {
        throw std::runtime_error(
            "Invalid critical pressure for material: " +
            material.name
            );
    }

    if (material.criticalVolumeM3PerKmol <= 0.0) {
        throw std::runtime_error(
            "Invalid critical volume for material: " +
            material.name
            );
    }

    if (material.molarMassKgPerKmol <= 0.0) {
        throw std::runtime_error(
            "Invalid molar mass for material: " +
            material.name
            );
    }
}

} // namespace

Material makeMaterialFromChemsep(
    const Component component,
    const ChemsepComponent& chemsepComponent
    )
{
    Material material;

    material.component =
        component;

    // Пока расчетная модель еще использует старые имена C2H6/C5H12/H2O.
    // Поэтому оставляем имя из enum, а не "Ethane"/"N-pentane"/"Water".
    material.name =
        componentName(component);

    material.criticalTemperatureK =
        requireFiniteScalar(
            chemsepComponent,
            "Tc_K"
            );

    material.criticalPressurePa =
        requireFiniteScalar(
            chemsepComponent,
            "Pc_Pa"
            );

    material.criticalVolumeM3PerKmol =
        requireFiniteScalar(
            chemsepComponent,
            "Vc_m3_per_kmol"
            );

    material.criticalCompressibility =
        requireFiniteScalar(
            chemsepComponent,
            "Zc"
            );

    material.acentricFactor =
        requireFiniteScalar(
            chemsepComponent,
            "omega"
            );

    material.molarMassKgPerKmol =
        requireFiniteScalar(
            chemsepComponent,
            "MW_kg_per_kmol"
            );

    material.antoine =
        makeVapourPressureCoefficients(
            chemsepComponent
            );

    validateMaterial(material);

    return material;
}

MaterialList makeMaterialListFromChemsep(
    const std::vector<Component>& components,
    const std::vector<const ChemsepComponent*>& chemsepComponents
    )
{
    if (components.size() != chemsepComponents.size()) {
        throw std::runtime_error(
            "Cannot build MaterialList: components size differs from ChemSep components size"
            );
    }

    if (components.empty()) {
        throw std::runtime_error(
            "Cannot build MaterialList: empty component list"
            );
    }

    MaterialList materials;
    materials.reserve(components.size());

    for (std::size_t i = 0; i < components.size(); ++i) {
        if (chemsepComponents[i] == nullptr) {
            throw std::runtime_error(
                "Cannot build MaterialList: null ChemSep component pointer"
                );
        }

        materials.push_back(
            makeMaterialFromChemsep(
                components[i],
                *chemsepComponents[i]
                )
            );
    }

    return materials;
}
