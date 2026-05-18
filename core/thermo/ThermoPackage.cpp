#include "ThermoPackage.h"

ThermoPackage::ThermoPackage()
    : materials_(createDefaultMaterials()),
    eos_(materials_),
    enthalpy_(materials_)
{
}

ThermoPackage::ThermoPackage(
    const MaterialList& materials
    )
    : materials_(materials),
    eos_(materials_),
    enthalpy_(materials_)
{
}

const EnthalpyModel& ThermoPackage::enthalpy() const
{
    return enthalpy_;
}

const MaterialList& ThermoPackage::materials() const
{
    return materials_;
}

const Material& ThermoPackage::material(
    Component component
    ) const
{
    return materials_[componentIndex(component)];
}

const PengRobinsonEOS& ThermoPackage::eos() const
{
    return eos_;
}

FlashResult ThermoPackage::flash(
    double pressurePa,
    double temperatureK,
    const Composition& zOverall,
    double tolerance,
    int maxIterations
    ) const
{
    FlashCalculation flashCalculation(
        materials_,
        eos_
        );

    FlashResult result =
        flashCalculation.calculate(
            pressurePa,
            temperatureK,
            zOverall,
            tolerance,
            maxIterations
            );

    result.molarEnthalpyLiquidJPerKmol =
        enthalpy_.phaseMolarEnthalpyJPerKmol(
            pressurePa,
            temperatureK,
            result.xLiquid,
            result.zLiquid
            );

    result.molarEnthalpyVaporJPerKmol =
        enthalpy_.phaseMolarEnthalpyJPerKmol(
            pressurePa,
            temperatureK,
            result.yVapor,
            result.zVapor
            );

    result.molarEnthalpyMixtureJPerKmol =
        (1.0 - result.beta)
            * result.molarEnthalpyLiquidJPerKmol
        + result.beta
              * result.molarEnthalpyVaporJPerKmol;

    return result;
}
