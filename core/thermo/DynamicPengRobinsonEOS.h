#pragma once

#include "../BinaryInteractionDatabase.h"
#include "../Component.h"

#include <cstddef>
#include <vector>

using DynamicComposition = std::vector<double>;

struct DynamicMixtureParameters {
    double a = 0.0;
    double b = 0.0;
    double A = 0.0;
    double B = 0.0;
};

struct DynamicZFactors {
    double liquid = 0.0;
    double vapor = 0.0;
    std::vector<double> roots;
};

class DynamicPengRobinsonEOS {
public:
    DynamicPengRobinsonEOS(
        ComponentSet components,
        BinaryInteractionMatrix kij
        );

    std::size_t componentCount() const;

    DynamicMixtureParameters computeMixtureParameters(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& composition
        ) const;

    DynamicZFactors computeZFactors(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& composition
        ) const;

    DynamicComposition computeFugacityCoefficients(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& composition,
        double zFactor
        ) const;

    double computeDepartureEnthalpyJPerKmol(
        double pressurePa,
        double temperatureK,
        const DynamicComposition& composition,
        double zFactor
        ) const;

    const ComponentProperties& component(
        std::size_t componentIndex
        ) const;

    double computePressureFromMolarVolume(
        double temperatureK,
        const DynamicComposition& composition,
        double molarVolumeM3PerKmol
        ) const;

    double computeZFactorFromMolarVolume(
        double pressurePa,
        double temperatureK,
        double molarVolumeM3PerKmol
        ) const;

private:
    ComponentSet components_;
    BinaryInteractionMatrix kij_;

    std::vector<double> bPure_;
    std::vector<double> mPure_;

    void computePureParameters();

    double computeAlpha(
        std::size_t componentIndex,
        double temperatureK
        ) const;

    std::vector<double> computeAPure(
        double temperatureK
        ) const;

    std::vector<double> solveCubicEquation(
        double c2,
        double c1,
        double c0
        ) const;

    void validateComposition(
        const DynamicComposition& composition
        ) const;
};
