#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

enum class Component : std::size_t {
    C2H6 = 0,
    C5H12 = 1,
    H2O = 2
};

constexpr std::size_t ComponentCount = 3;

using Composition = std::vector<double>;

inline Composition makeComposition(
    std::size_t size = ComponentCount,
    double value = 0.0
    )
{
    return Composition(size, value);
}

inline Composition makeLegacyComposition(
    double c2h6,
    double c5h12,
    double h2o
    )
{
    Composition composition(ComponentCount, 0.0);

    composition[0] = c2h6;
    composition[1] = c5h12;
    composition[2] = h2o;

    return composition;
}

inline void requireCompositionSize(
    const Composition& composition,
    std::size_t expectedSize,
    const std::string& context
    )
{
    if (composition.size() != expectedSize) {
        throw std::runtime_error(
            context +
            ": invalid composition size. Expected " +
            std::to_string(expectedSize) +
            ", got " +
            std::to_string(composition.size())
            );
    }
}

inline std::size_t componentIndex(Component component)
{
    return static_cast<std::size_t>(component);
}

inline std::string componentName(Component component)
{
    switch (component) {
    case Component::C2H6:
        return "C2H6";
    case Component::C5H12:
        return "C5H12";
    case Component::H2O:
        return "H2O";
    }

    return "Unknown";
}

inline double molecularWeightKgPerKmol(Component component)
{
    switch (component) {
    case Component::C2H6:
        return 30.06904;
    case Component::C5H12:
        return 72.14878;
    case Component::H2O:
        return 18.01528;
    }

    return 0.0;
}
