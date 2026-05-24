#include "ChemsepEnthalpyBuilder.h"

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

void validateCpRppCorrelation(
    const ChemsepComponent& component,
    const ChemsepCorrelation& cp
    )
{
    if (cp.eqno != 100) {
        std::ostringstream message;
        message
            << "Unsupported ChemSep ideal_gas_cp_rpp_eqno for "
            << component.componentName()
            << ". Expected 100, got "
            << cp.eqno;

        throw std::runtime_error(message.str());
    }

    if (
        !std::isfinite(cp.A) ||
        !std::isfinite(cp.B) ||
        !std::isfinite(cp.C) ||
        !std::isfinite(cp.D) ||
        !std::isfinite(cp.E)
        ) {
        throw std::runtime_error(
            "ChemSep ideal_gas_cp_rpp coefficients are incomplete for " +
            component.componentName()
            );
    }


}

void validateEnthalpyData(
    const ComponentEnthalpyData& data
    )
{
    if (!std::isfinite(data.referenceTemperatureK) || data.referenceTemperatureK <= 0.0) {
        throw std::runtime_error(
            "Invalid reference temperature for enthalpy data: " +
            data.name
            );
    }

    if (!std::isfinite(data.referenceMolarEnthalpyJPerKmol)) {
        throw std::runtime_error(
            "Invalid reference enthalpy for enthalpy data: " +
            data.name
            );
    }

    if (
        !std::isfinite(data.cpA) ||
        !std::isfinite(data.cpB) ||
        !std::isfinite(data.cpC) ||
        !std::isfinite(data.cpD) ||
        !std::isfinite(data.cpE)
        ) {
        throw std::runtime_error(
            "Invalid Cp coefficients for enthalpy data: " +
            data.name
            );
    }
}

} // namespace

ComponentEnthalpyData makeEnthalpyDataFromChemsep(
    const std::string& componentKey,
    const ChemsepComponent& chemsepComponent
    )
{
    const ChemsepCorrelation& cp =
        requireCorrelation(
            chemsepComponent,
            "ideal_gas_cp_rpp"
            );

    validateCpRppCorrelation(
        chemsepComponent,
        cp
        );

    ComponentEnthalpyData data;

    // Пока расчетная модель использует старые имена C2H6/C5H12/H2O.
    data.name =
        componentKey;

    data.referenceTemperatureK =
        298.15;

    data.referenceMolarEnthalpyJPerKmol =
        requireFiniteScalar(
            chemsepComponent,
            "h_form_ig_J_per_kmol"
            );

    // ideal_gas_cp_rpp_eqno = 100:
    // Cp[J/(kmol*K)] = A + B*T + C*T^2 + D*T^3 + E*T^4
    //
    // Текущий EnthalpyModel использует ту же форму:
    //
    // cpA = A
    // cpB = B
    // cpC = C
    // cpD = D
    // cpE = E

    data.cpA = cp.A;
    data.cpB = cp.B;
    data.cpC = cp.C;
    data.cpD = cp.D;
    data.cpE = cp.E;

    validateEnthalpyData(data);

    return data;
}

ComponentEnthalpyDataList makeEnthalpyDataListFromChemsep(
    const std::vector<std::string>& componentKeys,
    const std::vector<const ChemsepComponent*>& chemsepComponents
    )
{
    if (componentKeys.size() != chemsepComponents.size()) {
        throw std::runtime_error(
            "Cannot build ComponentEnthalpyDataList: componentKeys size differs from ChemSep components size"
            );
    }

    if (componentKeys.empty()) {
        throw std::runtime_error(
            "Cannot build ComponentEnthalpyDataList: empty component list"
            );
    }

    ComponentEnthalpyDataList enthalpyData;
    enthalpyData.reserve(componentKeys.size());

    for (std::size_t i = 0; i < componentKeys.size(); ++i) {
        if (chemsepComponents[i] == nullptr) {
            throw std::runtime_error(
                "Cannot build ComponentEnthalpyDataList: null ChemSep component pointer"
                );
        }

        enthalpyData.push_back(
            makeEnthalpyDataFromChemsep(
                componentKeys[i],
                *chemsepComponents[i]
                )
            );
    }

    return enthalpyData;
}
