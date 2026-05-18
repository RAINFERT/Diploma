#include "ReactorState.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>

ReactorState::ReactorState()
    : massesKg_{0.0, 1.0, 500.0},
    energyJ_(0.0),
    temperatureC_(100.0),
    pressureBar_(2.5)
{
}

double ReactorState::massKg(Component component) const
{
    return massesKg_[componentIndex(component)];
}

void ReactorState::setMassKg(Component component, double value)
{
    if (value < 0.0)
    {
        throw std::invalid_argument("Mass cannot be negative");
    }

    massesKg_[componentIndex(component)] = value;
}

void ReactorState::addMassKg(
    Component component,
    double deltaKg
    )
{
    const double newValue =
        massKg(component) + deltaKg;

    if (newValue < 0.0)
    {
        throw std::runtime_error("Mass became negative");
    }

    massesKg_[componentIndex(component)] = newValue;
}

double ReactorState::molesKmol(Component component) const
{
    const double mass = massKg(component);
    const double molecularWeight = molecularWeightKgPerKmol(component);

    return mass / molecularWeight;
}

double ReactorState::totalMolesKmol() const
{
    double total = 0.0;

    total += molesKmol(Component::C2H6);
    total += molesKmol(Component::C5H12);
    total += molesKmol(Component::H2O);

    return total;
}

double ReactorState::moleFraction(Component component) const
{
    const double total = totalMolesKmol();

    if (total <= 0.0)
    {
        return 0.0;
    }

    return molesKmol(component) / total;
}

std::array<double, ComponentCount> ReactorState::moleFractions() const
{
    return {
        moleFraction(Component::C2H6),
        moleFraction(Component::C5H12),
        moleFraction(Component::H2O)
    };
}

double ReactorState::energyJ() const
{
    return energyJ_;
}

void ReactorState::setEnergyJ(double value)
{
    energyJ_ = value;
}

double ReactorState::temperatureC() const
{
    return temperatureC_;
}

void ReactorState::setTemperatureC(double value)
{
    temperatureC_ = value;
}

double ReactorState::pressureBar() const
{
    return pressureBar_;
}

void ReactorState::setPressureBar(double value)
{
    if (value <= 0.0)
    {
        throw std::invalid_argument("Pressure must be positive");
    }

    pressureBar_ = value;
}

std::string ReactorState::toString() const
{
    std::ostringstream out;

    out << std::fixed << std::setprecision(6);

    out << "Reactor state:\n";

    out << "  Masses:\n";
    out << "    M(C2H6)  = " << massKg(Component::C2H6) << " kg\n";
    out << "    M(C5H12) = " << massKg(Component::C5H12) << " kg\n";
    out << "    M(H2O)   = " << massKg(Component::H2O) << " kg\n";

    out << "  Amounts:\n";
    out << "    n(C2H6)  = " << molesKmol(Component::C2H6) << " kmol\n";
    out << "    n(C5H12) = " << molesKmol(Component::C5H12) << " kmol\n";
    out << "    n(H2O)   = " << molesKmol(Component::H2O) << " kmol\n";
    out << "    n(total) = " << totalMolesKmol() << " kmol\n";

    out << "  Mole fractions:\n";
    out << "    Z(C2H6)  = " << moleFraction(Component::C2H6) << "\n";
    out << "    Z(C5H12) = " << moleFraction(Component::C5H12) << "\n";
    out << "    Z(H2O)   = " << moleFraction(Component::H2O) << "\n";

    out << "  Other variables:\n";
    out << "    E = " << energyJ_ << " J\n";
    out << "    T = " << temperatureC_ << " C\n";
    out << "    P = " << pressureBar_ << " bar\n";

    return out.str();
}
