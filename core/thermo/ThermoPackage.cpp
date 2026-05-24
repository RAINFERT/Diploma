#include "ThermoPackage.h"
#include "GibbsFlashCalculation.h"
#include <iostream>

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
    GibbsFlashCalculation  flashCalculation(
        materials_,
        eos_
        );

    FlashResult result =
        flashCalculation.calculate(
            pressurePa,
            temperatureK,
            zOverall,
            tolerance,
            std::max(maxIterations, 80)
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

    static int flashCallCounter = 0;
    ++flashCallCounter;

    /*
    if (flashCallCounter <= 20 || flashCallCounter % 100 == 0)
    {
        std::cout
            << "[flash #" << flashCallCounter << "] "
            << "method = " << flashMethodToString(result.method)
            << ", status = " << flashStatusToString(result.status)
            << ", beta = " << result.beta
            << ", accepted = " << result.diagnostics.acceptedRadauSteps
            << ", rejected = " << result.diagnostics.rejectedRadauSteps
            << "\n";
    }

    */

    if (result.status == FlashStatus::NotConverged)
    {
        std::cerr
            << "[WARNING] Flash did not converge. "
            << "P = " << pressurePa
            << ", T = " << temperatureK
            << ", beta = " << result.beta
            << ", method = " << flashMethodToString(result.method)
            << "\n";
    }

    return result;
}
