#pragma once

#include "Component.h"

#include <array>
#include <string>

class ReactorState
{
public:
    ReactorState();

    double massKg(Component component) const;
    void setMassKg(Component component, double value);
    void addMassKg(
        Component component,
        double deltaKg
        );

    double molesKmol(Component component) const;
    double totalMolesKmol() const;

    double moleFraction(Component component) const;
    Composition moleFractions() const;

    double energyJ() const;
    void setEnergyJ(double value);

    double temperatureC() const;
    void setTemperatureC(double value);

    double pressureBar() const;
    void setPressureBar(double value);

    std::string toString() const;

private:
    Composition massesKg_;

    double energyJ_;
    double temperatureC_;
    double pressureBar_;
};
