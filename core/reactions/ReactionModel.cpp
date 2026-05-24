#include "ReactionModel.h"

#include <stdexcept>

ReactionModel::ReactionModel(
    const MaterialList& materials
    )
    : materials_(materials)
{
}

void ReactionModel::setForwardRateConstant(
    double value
    )
{
    if (value < 0.0)
    {
        throw std::invalid_argument("Forward reaction rate constant cannot be negative");
    }

    kForwardPerS_ = value;
}

void ReactionModel::setReverseRateConstant(
    double value
    )
{
    if (value < 0.0)
    {
        throw std::invalid_argument("Reverse reaction rate constant cannot be negative");
    }

    kReversePerS_ = value;
}

void ReactionModel::setEnabled(const bool enabled)
{
    enabled_ = enabled;
}

bool ReactionModel::enabled() const
{
    return enabled_;
}

double ReactionModel::forwardRateConstant() const
{
    return kForwardPerS_;
}

double ReactionModel::reverseRateConstant() const
{
    return kReversePerS_;
}

ReactionRates ReactionModel::computeRates(
    double temperatureK,
    double pressurePa,
    const Composition& liquidComposition,
    double liquidMolarDensityKmolPerM3
    ) const
{
    requireCompositionSize(
        liquidComposition,
        materials_.size(),
        "ReactionModel::computeRates liquidComposition"
        );

    ReactionRates result;
    result.componentRatesKmolPerM3S =
        makeComposition(materials_.size());

    if (!enabled_) {
        return result;
    }

    throw std::runtime_error(
        "ReactionModel is not implemented for arbitrary ChemSep component sets. "
        "Set reactions.enabled=false or implement stoichiometric reactions in config."
        );

    if (temperatureK <= 0.0)
    {
        throw std::invalid_argument("Temperature must be positive");
    }

    if (pressurePa <= 0.0)
    {
        throw std::invalid_argument("Pressure must be positive");
    }

    if (liquidMolarDensityKmolPerM3 < 0.0)
    {
        throw std::invalid_argument("Liquid molar density cannot be negative");
    }

    if (materials_.size() != ComponentCount) {
        throw std::runtime_error(
            "Legacy reaction model supports only C2H6/C5H12/H2O. Disable reactions for arbitrary component sets."
            );
    }

    const std::size_t c2Index =
        componentIndex(Component::C2H6);

    const std::size_t c5Index =
        componentIndex(Component::C5H12);

    const double concentrationC2 =
        liquidComposition[c2Index]
        * liquidMolarDensityKmolPerM3;

    const double concentrationC5 =
        liquidComposition[c5Index]
        * liquidMolarDensityKmolPerM3;

    const double rForward =
        kForwardPerS_
        * concentrationC2;

    const double rReverse =
        kReversePerS_
        * concentrationC5;

    result.forwardRateKmolPerM3S = rForward;
    result.reverseRateKmolPerM3S = rReverse;

    const double mwC2 =
        materials_[c2Index].molarMassKgPerKmol;

    const double mwC5 =
        materials_[c5Index].molarMassKgPerKmol;

    const double massConservingRatio =
        mwC5 / mwC2;

    result.componentRatesKmolPerM3S[c2Index] =
        massConservingRatio
        * (-rForward + rReverse);

    result.componentRatesKmolPerM3S[c5Index] =
        rForward - rReverse;

    result.componentRatesKmolPerM3S[componentIndex(Component::H2O)] =
        0.0;

    return result;
}
