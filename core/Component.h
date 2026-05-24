#pragma once

#include <cstddef>
#include <string>

enum class Component : std::size_t
{
    C2H6 = 0,
    C5H12 = 1,
    H2O = 2
};

constexpr std::size_t ComponentCount = 3;

using Composition = std::array<double, ComponentCount>;

inline std::size_t componentIndex(Component component)
{
    return static_cast<std::size_t>(component);
}

inline std::string componentName(Component component)
{
    switch (component)
    {
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
    switch (component)
    {
    case Component::C2H6:
        return 30.06904;

    case Component::C5H12:
        return 72.14878;

    case Component::H2O:
        return 18.01528;
    }

    return 0.0;
}
