#include "ReactorState.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>

ReactorState::ReactorState()
    : ReactorState(
          std::vector<std::string>{"C2H6", "C5H12", "H2O"},
          std::vector<double>{
              molecularWeightKgPerKmol(Component::C2H6),
              molecularWeightKgPerKmol(Component::C5H12),
              molecularWeightKgPerKmol(Component::H2O)
          }
          )
{
    setMassKg(componentIndex(Component::C2H6), 0.0);
    setMassKg(componentIndex(Component::C5H12), 1.0);
    setMassKg(componentIndex(Component::H2O), 500.0);
}

ReactorState::ReactorState(
    std::vector<std::string> componentNames,
    std::vector<double> molarMassesKgPerKmol
    )
    : massesKg_(makeComposition(componentNames.size())),
    energyJ_(0.0),
    temperatureC_(100.0),
    pressureBar_(2.5),
    componentNames_(std::move(componentNames)),
    molarMassesKgPerKmol_(std::move(molarMassesKgPerKmol))
{
    if (componentNames_.empty()) {
        throw std::runtime_error(
            "ReactorState must contain at least one component"
            );
    }

    if (componentNames_.size() != molarMassesKgPerKmol_.size()) {
        throw std::runtime_error(
            "ReactorState component names size differs from molar masses size"
            );
    }
}

std::size_t ReactorState::componentCount() const
{
    return massesKg_.size();
}

const std::string& ReactorState::componentName(
    const std::size_t index
    ) const
{
    return componentNames_.at(index);
}

double ReactorState::molarMassKgPerKmol(
    const std::size_t index
    ) const
{
    return molarMassesKgPerKmol_.at(index);
}

double ReactorState::massKg(
    const std::size_t index
    ) const
{
    return massesKg_.at(index);
}

void ReactorState::setMassKg(
    const std::size_t index,
    const double value
    )
{
    if (value < 0.0) {
        throw std::invalid_argument(
            "Component mass cannot be negative"
            );
    }

    massesKg_.at(index) = value;
}

void ReactorState::addMassKg(
    const std::size_t index,
    const double deltaKg
    )
{
    const double updated =
        massKg(index) + deltaKg;

    if (updated < 0.0) {
        throw std::runtime_error(
            "Component mass became negative"
            );
    }

    setMassKg(index, updated);
}

double ReactorState::molesKmol(
    const std::size_t index
    ) const
{
    const double molarMass =
        molarMassKgPerKmol(index);

    if (molarMass <= 0.0) {
        throw std::runtime_error(
            "Molar mass must be positive"
            );
    }

    return massKg(index) / molarMass;
}

double ReactorState::moleFraction(
    const std::size_t index
    ) const
{
    const double total =
        totalMolesKmol();

    if (total <= 0.0) {
        return 0.0;
    }

    return molesKmol(index) / total;
}

double ReactorState::massKg(Component component) const
{
    return massKg(componentIndex(component));
}

void ReactorState::setMassKg(Component component, double value)
{
    setMassKg(componentIndex(component), value);
}

void ReactorState::addMassKg(Component component, double deltaKg)
{
    addMassKg(componentIndex(component), deltaKg);
}

double ReactorState::molesKmol(Component component) const
{
    return molesKmol(componentIndex(component));
}

double ReactorState::moleFraction(Component component) const
{
    return moleFraction(componentIndex(component));
}

double ReactorState::totalMolesKmol() const
{
    double total = 0.0;

    for (std::size_t i = 0; i < componentCount(); ++i) {
        total += molesKmol(i);
    }

    return total;
}

Composition ReactorState::moleFractions() const
{
    Composition fractions =
        makeComposition(componentCount());

    const double total =
        totalMolesKmol();

    if (total <= 0.0) {
        return fractions;
    }

    for (std::size_t i = 0; i < componentCount(); ++i) {
        fractions[i] =
            molesKmol(i) / total;
    }

    return fractions;
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
    for (std::size_t i = 0; i < componentCount(); ++i) {
        out << "    M(" << componentName(i) << ") = "
            << massesKg_[i] << " kg\n";
    }

    out << "  Amounts:\n";
    for (std::size_t i = 0; i < componentCount(); ++i) {
        out << "    n(" << componentName(i) << ") = "
            << molesKmol[i] << " kmol\n";
    }
    out << "    n(total) = " << totalMolesKmol() << " kmol\n";

    out << "  Mole fractions:\n";
    for (std::size_t i = 0; i < componentCount(); ++i) {
        out << "    n(" << componentName(i) << ") = "
            << moleFraction[i] << "\n";
    }

    out << "  Other variables:\n";
    out << "    E = " << energyJ_ << " J\n";
    out << "    T = " << temperatureC_ << " C\n";
    out << "    P = " << pressureBar_ << " bar\n";

    return out.str();
}
