#pragma once

#include "Material.h"

#include <array>
#include <cstddef>

using Composition = std::array<double, ComponentCount>;

struct MixtureParameters
{
    double a;
    double b;
};

struct ZFactors
{
    double liquid;
    double middle;
    double vapor;

    std::size_t rootCount;
};

class PengRobinsonEOS
{
public:
    explicit PengRobinsonEOS(const MaterialList& materials);

    MixtureParameters computeMixtureParameters(
        const Composition& composition,
        double temperatureK
        ) const;

    ZFactors computeZFactors(
        double pressurePa,
        double temperatureK,
        const Composition& composition
        ) const;

    Composition computeFugacityCoefficients(
        double pressurePa,
        double temperatureK,
        const Composition& composition,
        double zFactor
        ) const;

private:
    static constexpr double R = 8.314462618;

    MaterialList materials_;

    std::array<double, ComponentCount> bi_;
    std::array<double, ComponentCount> mi_;

    std::array<std::array<double, ComponentCount>, ComponentCount> kij_;

    void computePureParameters();

    double computeAlpha(
        std::size_t componentIndex,
        double temperatureK
        ) const;

    std::array<double, ComponentCount> computeAi(
        double temperatureK
        ) const;
};
