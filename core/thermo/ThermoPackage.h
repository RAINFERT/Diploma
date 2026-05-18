#pragma once

#include "FlashCalculation.h"
#include "Material.h"
#include "PengRobinsonEOS.h"
#include "EnthalpyModel.h"

class ThermoPackage
{
public:
    ThermoPackage();

    explicit ThermoPackage(
        const MaterialList& materials
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
};
