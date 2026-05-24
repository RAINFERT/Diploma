#pragma once

#include "Material.h"

#include <array>
#include <cstddef>

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

    double pressureFromMolarVolume(
        double temperatureK,
        double molarVolumeM3PerMol,
        const Composition& composition
        ) const;

    Composition computeLogFugacityCoefficients(
        double pressurePa,
        double temperatureK,
        const Composition& composition,
        double zFactor
        ) const;

    double pseudoCriticalMolarVolume(
        const Composition& composition
        ) const;

    double pseudoCriticalMolarVolumeM3PerMol(
        const Composition& composition
        ) const;



private:
    static constexpr double R = 8.314462618;

    MaterialList materials_;

    std::vector<double> bi_;
    std::vector<double> mi_;
    std::vector<std::vector<double>> kij_;

    void computePureParameters();

    double computeAlpha(
        std::size_t componentIndex,
        double temperatureK
        ) const;

    std::vector<double> computeAi(double temperatureK) const;
};
