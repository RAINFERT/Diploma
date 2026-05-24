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

    static int totalFlashCalls = 0;
    static int rrBoundaryCalls = 0;
    static int rrFullPrecheckCalls = 0;
    static int gibbsRadauCalls = 0;
    static int singlePrecheckCalls = 0;
    static int otherCalls = 0;

    ++totalFlashCalls;

    switch (result.method)
    {
    case FlashMethod::HybridBoundaryRachfordRice:
        ++rrBoundaryCalls;
        break;

    case FlashMethod::HybridFullRachfordRicePrecheck:
        ++rrFullPrecheckCalls;
        break;

    case FlashMethod::GibbsRadau:
        ++gibbsRadauCalls;
        break;

    case FlashMethod::HybridSinglePhasePrecheck:
        ++singlePrecheckCalls;
        break;

    default:
        ++otherCalls;
        break;
    }

    if (totalFlashCalls % 1000 == 0)
    {
        std::cout
            << "[flash stats] total=" << totalFlashCalls
            << ", boundaryRR=" << rrBoundaryCalls
            << ", fullRRprecheck=" << rrFullPrecheckCalls
            << ", GibbsRadau=" << gibbsRadauCalls
            << ", singlePrecheck=" << singlePrecheckCalls
            << ", other=" << otherCalls
            << "\n";
    }

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
