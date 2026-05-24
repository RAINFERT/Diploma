#pragma once

#include "FlashCalculation.h"
#include "Material.h"
#include "PengRobinsonEOS.h"
#include "EnthalpyModel.h"

enum class ThermoFlashMode
{
    RachfordRice,
    GibbsRadauHybrid
};

inline const char* thermoFlashModeToString(ThermoFlashMode mode)
{
    switch (mode)
    {
    case ThermoFlashMode::RachfordRice:
        return "Rachford-Rice";

    case ThermoFlashMode::GibbsRadauHybrid:
        return "Gibbs-Radau hybrid";

    default:
        return "Unknown";
    }
}

class ThermoPackage
{
public:
    ThermoPackage();

    explicit ThermoPackage(
        ThermoFlashMode flashMode
        );

    void setFlashMode(
        ThermoFlashMode flashMode
        );

    ThermoFlashMode flashMode() const;

    explicit ThermoPackage(
        const MaterialList& materials,
        ThermoFlashMode flashMode = ThermoFlashMode::GibbsRadauHybrid
        );

    ThermoPackage(
        const MaterialList& materials,
        const ComponentEnthalpyDataList& enthalpyData,
        ThermoFlashMode flashMode = ThermoFlashMode::GibbsRadauHybrid
        );

    const EnthalpyModel& enthalpy() const;

    const MaterialList& materials() const;

    const Material& material(
        Component component
        ) const;

    const PengRobinsonEOS& eos() const;

    FlashResult flash(
        double pressurePa,
        double temperatureK,
        const Composition& zOverall,
        double tolerance = 1.0e-8,
        int maxIterations = 100
        ) const;

private:
    MaterialList materials_;
    PengRobinsonEOS eos_;
    EnthalpyModel enthalpy_;

    ThermoFlashMode flashMode_ = ThermoFlashMode::GibbsRadauHybrid;
};
