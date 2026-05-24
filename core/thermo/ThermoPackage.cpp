#include "ThermoPackage.h"
#include "FlashCalculation.h"
#include "GibbsFlashCalculation.h"

#include <algorithm>
#include <iostream>

namespace
{
constexpr bool PrintFlashWarnings = true;
}

ThermoPackage::ThermoPackage()
    : ThermoPackage(
        createDefaultMaterials(),
        ThermoFlashMode::GibbsRadauHybrid
    )
{
}

ThermoPackage::ThermoPackage(
    ThermoFlashMode flashMode
)
    : ThermoPackage(
        createDefaultMaterials(),
        flashMode
    )
{
}

ThermoPackage::ThermoPackage(
    const MaterialList& materials,
    ThermoFlashMode flashMode
)
    : materials_(materials),
      eos_(materials_),
      enthalpy_(materials_),
      flashMode_(flashMode)
{
}

void ThermoPackage::setFlashMode(
    ThermoFlashMode flashMode
)
{
    flashMode_ = flashMode;
}

ThermoFlashMode ThermoPackage::flashMode() const
{
    return flashMode_;
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
    FlashResult result;

    if (flashMode_ == ThermoFlashMode::RachfordRice)
    {
        FlashCalculation flashCalculation(
            materials_,
            eos_
        );

        result =
            flashCalculation.calculate(
                pressurePa,
                temperatureK,
                zOverall,
                tolerance,
                maxIterations
            );

        result.method = FlashMethod::RachfordRice;
    }
    else
    {
        GibbsFlashCalculation flashCalculation(
            materials_,
            eos_
        );

        result =
            flashCalculation.calculate(
                pressurePa,
                temperatureK,
                zOverall,
                tolerance,
                std::max(maxIterations, 80)
            );
    }

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

    if (
        PrintFlashWarnings
        &&
        result.status == FlashStatus::NotConverged
    )
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
