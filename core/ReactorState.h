#pragma once

#include "Component.h"

#include <array>
#include <string>

class ReactorState
{
public:
    ReactorState();

    ReactorState(
        std::vector<std::string> componentNames,
        std::vector<double> molarMassesKgPerKmol
        );

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

    std::size_t componentCount() const;

    const std::string& componentName(
        std::size_t index
        ) const;

    double molarMassKgPerKmol(
        std::size_t index
        ) const;

    double massKg(
        std::size_t index
        ) const;

    void setMassKg(
        std::size_t index,
        double value
        );

    void addMassKg(
        std::size_t index,
        double deltaKg
        );

    double molesKmol(
        std::size_t index
        ) const;

    double moleFraction(
        std::size_t index
        ) const;

private:
    Composition massesKg_;

    std::vector<std::string> componentNames_;
    std::vector<double> molarMassesKgPerKmol_;

    double energyJ_;
    double temperatureC_;
    double pressureBar_;
};
