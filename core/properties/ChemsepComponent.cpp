#include "ChemsepComponent.h"

#include <cmath>

namespace {

const std::string emptyString;

} // namespace

void ChemsepComponent::setField(
    const std::string& name,
    const std::string& value
    )
{
    rawFields_[name] = value;
}

void ChemsepComponent::setScalar(
    const std::string& name,
    const double value
    )
{
    scalarFields_[name] = value;
}

void ChemsepComponent::setCorrelation(
    const std::string& name,
    const ChemsepCorrelation& correlation
    )
{
    correlations_[name] = correlation;
}

const std::string& ChemsepComponent::field(
    const std::string& name
    ) const
{
    const auto it = rawFields_.find(name);

    if (it == rawFields_.end()) {
        return emptyString;
    }

    return it->second;
}

bool ChemsepComponent::hasField(
    const std::string& name
    ) const
{
    return rawFields_.find(name) != rawFields_.end();
}

bool ChemsepComponent::hasScalar(
    const std::string& name
    ) const
{
    return scalarFields_.find(name) != scalarFields_.end();
}

double ChemsepComponent::scalar(
    const std::string& name
    ) const
{
    const auto it = scalarFields_.find(name);

    if (it == scalarFields_.end()) {
        throw std::runtime_error(
            "ChemSep scalar property is missing: " + name
            );
    }

    return it->second;
}

double ChemsepComponent::scalarOrNaN(
    const std::string& name
    ) const
{
    const auto it = scalarFields_.find(name);

    if (it == scalarFields_.end()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return it->second;
}

bool ChemsepComponent::hasCorrelation(
    const std::string& name
    ) const
{
    return correlations_.find(name) != correlations_.end();
}

const ChemsepCorrelation& ChemsepComponent::correlation(
    const std::string& name
    ) const
{
    const auto it = correlations_.find(name);

    if (it == correlations_.end()) {
        throw std::runtime_error(
            "ChemSep correlation is missing: " + name
            );
    }

    return it->second;
}

const std::unordered_map<std::string, std::string>&
ChemsepComponent::rawFields() const
{
    return rawFields_;
}

const std::unordered_map<std::string, double>&
ChemsepComponent::scalarFields() const
{
    return scalarFields_;
}

const std::unordered_map<std::string, ChemsepCorrelation>&
ChemsepComponent::correlations() const
{
    return correlations_;
}

std::string ChemsepComponent::componentName() const
{
    return field("component_name");
}

std::string ChemsepComponent::formulaChemsep() const
{
    return field("formula_chemsep");
}

std::string ChemsepComponent::casNumber() const
{
    return field("cas_number");
}

double ChemsepComponent::molarMassKgPerKmol() const
{
    return scalar("MW_kg_per_kmol");
}

double ChemsepComponent::criticalTemperatureK() const
{
    return scalar("Tc_K");
}

double ChemsepComponent::criticalPressurePa() const
{
    return scalar("Pc_Pa");
}

double ChemsepComponent::criticalVolumeM3PerKmol() const
{
    return scalar("Vc_m3_per_kmol");
}

double ChemsepComponent::criticalCompressibility() const
{
    return scalar("Zc");
}

double ChemsepComponent::acentricFactor() const
{
    return scalar("omega");
}

double ChemsepComponent::idealGasFormationEnthalpyJPerKmol() const
{
    return scalar("h_form_ig_J_per_kmol");
}
