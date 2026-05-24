#include "DynamicPengRobinsonEOS.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include <stdexcept>

namespace {

constexpr double GasConstantJPerKmolK = 8314.46261815324;
constexpr double OmegaA = 0.45724;
constexpr double OmegaB = 0.07780;
constexpr double Sqrt2 = 1.4142135623730951;

double clampValue(double value, double minValue, double maxValue) {
    return std::max(minValue, std::min(value, maxValue));
}

DynamicComposition normalizedComposition(
    const DynamicComposition& composition
    ) {
    double sum = std::accumulate(
        composition.begin(),
        composition.end(),
        0.0
        );

    if (sum <= 0.0) {
        throw std::runtime_error("Composition sum must be positive");
    }

    DynamicComposition result = composition;

    for (double& value : result) {
        if (value < 0.0) {
            throw std::runtime_error("Composition contains negative value");
        }

        value /= sum;
    }

    return result;
}

} // namespace

DynamicPengRobinsonEOS::DynamicPengRobinsonEOS(
    ComponentSet components,
    BinaryInteractionMatrix kij
    )
    : components_(std::move(components)),
    kij_(std::move(kij))
{
    const std::size_t n = components_.size();

    if (n == 0) {
        throw std::runtime_error("DynamicPengRobinsonEOS: empty component set");
    }

    if (kij_.size() != n) {
        throw std::runtime_error("DynamicPengRobinsonEOS: invalid kij row count");
    }

    for (const auto& row : kij_) {
        if (row.size() != n) {
            throw std::runtime_error("DynamicPengRobinsonEOS: invalid kij column count");
        }
    }

    bPure_.assign(n, 0.0);
    mPure_.assign(n, 0.0);

    computePureParameters();
}

std::size_t DynamicPengRobinsonEOS::componentCount() const {
    return components_.size();
}

void DynamicPengRobinsonEOS::validateComposition(
    const DynamicComposition& composition
    ) const {
    if (composition.size() != components_.size()) {
        throw std::runtime_error(
            "Composition size does not match component count"
            );
    }
}

void DynamicPengRobinsonEOS::computePureParameters() {
    const std::size_t n = components_.size();

    for (std::size_t i = 0; i < n; ++i) {
        const ComponentProperties& component = components_[i];

        const double tc = component.criticalTemperatureK;
        const double pc = component.criticalPressurePa;
        const double omega = component.acentricFactor;

        if (tc <= 0.0 || pc <= 0.0) {
            throw std::runtime_error(
                "Invalid critical properties for component: " +
                component.name
                );
        }

        bPure_[i] =
            OmegaB * GasConstantJPerKmolK * tc / pc;

        if (omega <= 0.491) {
            mPure_[i] =
                0.37464
                + 1.54226 * omega
                - 0.26992 * omega * omega;
        } else {
            mPure_[i] =
                0.379642
                + 1.48503 * omega
                - 0.164423 * omega * omega
                + 0.016666 * omega * omega * omega;
        }
    }
}

double DynamicPengRobinsonEOS::computeAlpha(
    std::size_t componentIndex,
    double temperatureK
    ) const {
    const ComponentProperties& component = components_[componentIndex];

    const double tr =
        temperatureK / component.criticalTemperatureK;

    if (tr <= 0.0) {
        throw std::runtime_error("Reduced temperature must be positive");
    }

    const double sqrtTr = std::sqrt(tr);

    const double factor =
        1.0 + mPure_[componentIndex] * (1.0 - sqrtTr);

    return factor * factor;
}

std::vector<double> DynamicPengRobinsonEOS::computeAPure(
    double temperatureK
    ) const {
    const std::size_t n = components_.size();

    std::vector<double> aPure(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        const ComponentProperties& component = components_[i];

        const double tc = component.criticalTemperatureK;
        const double pc = component.criticalPressurePa;
        const double alpha = computeAlpha(i, temperatureK);

        aPure[i] =
            OmegaA
            * GasConstantJPerKmolK
            * GasConstantJPerKmolK
            * tc
            * tc
            / pc
            * alpha;
    }

    return aPure;
}

DynamicMixtureParameters DynamicPengRobinsonEOS::computeMixtureParameters(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& composition
    ) const {
    validateComposition(composition);

    if (pressurePa <= 0.0 || temperatureK <= 0.0) {
        throw std::runtime_error("Pressure and temperature must be positive");
    }

    const DynamicComposition x = normalizedComposition(composition);
    const std::vector<double> aPure = computeAPure(temperatureK);

    const std::size_t n = components_.size();

    double aMix = 0.0;
    double bMix = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        bMix += x[i] * bPure_[i];
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const double aij =
                std::sqrt(aPure[i] * aPure[j])
                * (1.0 - kij_[i][j]);

            aMix += x[i] * x[j] * aij;
        }
    }

    DynamicMixtureParameters params;

    params.a = aMix;
    params.b = bMix;
    params.A =
        aMix * pressurePa
        / (GasConstantJPerKmolK * GasConstantJPerKmolK
           * temperatureK * temperatureK);

    params.B =
        bMix * pressurePa
        / (GasConstantJPerKmolK * temperatureK);

    return params;
}

std::vector<double> DynamicPengRobinsonEOS::solveCubicEquation(
    double c2,
    double c1,
    double c0
    ) const {
    // Решаем уравнение:
    // Z^3 + c2 Z^2 + c1 Z + c0 = 0

    const double pi = std::acos(-1.0);

    const double a = c2;
    const double b = c1;
    const double c = c0;

    const double p = b - a * a / 3.0;
    const double q =
        2.0 * a * a * a / 27.0
        - a * b / 3.0
        + c;

    const double discriminant =
        q * q / 4.0
        + p * p * p / 27.0;

    std::vector<double> roots;

    if (discriminant > 1.0e-14) {
        const double sqrtD = std::sqrt(discriminant);

        const double u = std::cbrt(-q / 2.0 + sqrtD);
        const double v = std::cbrt(-q / 2.0 - sqrtD);

        roots.push_back(u + v - a / 3.0);
    } else {
        if (std::abs(p) < 1.0e-14) {
            roots.push_back(-a / 3.0);
        } else {
            const double argument =
                clampValue(
                    (3.0 * q / (2.0 * p)) * std::sqrt(-3.0 / p),
                    -1.0,
                    1.0
                    );

            const double phi = std::acos(argument);
            const double radius = 2.0 * std::sqrt(-p / 3.0);

            for (int k = 0; k < 3; ++k) {
                const double root =
                    radius * std::cos((phi + 2.0 * pi * k) / 3.0)
                    - a / 3.0;

                roots.push_back(root);
            }
        }
    }

    roots.erase(
        std::remove_if(
            roots.begin(),
            roots.end(),
            [](double value) {
                return !std::isfinite(value) || value <= 0.0;
            }
            ),
        roots.end()
        );

    std::sort(roots.begin(), roots.end());

    roots.erase(
        std::unique(
            roots.begin(),
            roots.end(),
            [](double lhs, double rhs) {
                return std::abs(lhs - rhs) < 1.0e-10;
            }
            ),
        roots.end()
        );

    if (roots.empty()) {
        throw std::runtime_error("No positive real Z roots found");
    }

    return roots;
}

DynamicZFactors DynamicPengRobinsonEOS::computeZFactors(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& composition
    ) const {
    const DynamicMixtureParameters params =
        computeMixtureParameters(
            pressurePa,
            temperatureK,
            composition
            );

    const double A = params.A;
    const double B = params.B;

    // Peng-Robinson cubic:
    // Z^3 - (1 - B) Z^2 + (A - 3B^2 - 2B) Z
    // - (AB - B^2 - B^3) = 0

    const double c2 = -(1.0 - B);
    const double c1 = A - 3.0 * B * B - 2.0 * B;
    const double c0 = -(A * B - B * B - B * B * B);

    DynamicZFactors result;

    result.roots = solveCubicEquation(c2, c1, c0);

    result.liquid = result.roots.front();
    result.vapor = result.roots.back();

    return result;
}

DynamicComposition DynamicPengRobinsonEOS::computeFugacityCoefficients(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& composition,
    double zFactor
    ) const {
    validateComposition(composition);

    if (zFactor <= 0.0) {
        throw std::runtime_error("Z factor must be positive");
    }

    const DynamicComposition x = normalizedComposition(composition);

    const DynamicMixtureParameters params =
        computeMixtureParameters(
            pressurePa,
            temperatureK,
            x
            );

    const std::vector<double> aPure = computeAPure(temperatureK);

    const std::size_t n = components_.size();

    DynamicComposition phi(n, 1.0);

    const double aMix = params.a;
    const double bMix = params.b;
    const double A = params.A;
    const double B = params.B;

    if (aMix <= 0.0 || bMix <= 0.0 || B <= 1.0e-14) {
        return phi;
    }

    if (zFactor <= B) {
        throw std::runtime_error("Invalid Z factor: Z <= B");
    }

    const double logZMinusB = std::log(zFactor - B);

    const double logTermArgument =
        (zFactor + (1.0 + Sqrt2) * B)
        / (zFactor + (1.0 - Sqrt2) * B);

    if (logTermArgument <= 0.0) {
        throw std::runtime_error("Invalid logarithm argument in fugacity calculation");
    }

    const double logTerm = std::log(logTermArgument);

    for (std::size_t i = 0; i < n; ++i) {
        double sumAij = 0.0;

        for (std::size_t j = 0; j < n; ++j) {
            const double aij =
                std::sqrt(aPure[i] * aPure[j])
                * (1.0 - kij_[i][j]);

            sumAij += x[j] * aij;
        }

        const double biOverB = bPure_[i] / bMix;

        const double attractionTerm =
            2.0 * sumAij / aMix - biOverB;

        const double lnPhi =
            biOverB * (zFactor - 1.0)
            - logZMinusB
            - A / (2.0 * Sqrt2 * B)
                  * attractionTerm
                  * logTerm;

        phi[i] = std::exp(lnPhi);
    }

    return phi;
}

double DynamicPengRobinsonEOS::computeDepartureEnthalpyJPerKmol(
    double pressurePa,
    double temperatureK,
    const DynamicComposition& composition,
    double zFactor
    ) const {
    validateComposition(composition);

    if (pressurePa <= 0.0 || temperatureK <= 0.0) {
        throw std::runtime_error(
            "DynamicPengRobinsonEOS: pressure and temperature must be positive"
            );
    }

    if (zFactor <= 0.0) {
        throw std::runtime_error(
            "DynamicPengRobinsonEOS: Z factor must be positive"
            );
    }

    const DynamicComposition x =
        normalizedComposition(composition);

    const std::size_t n =
        components_.size();

    std::vector<double> aPure(n, 0.0);
    std::vector<double> daPureDt(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        const ComponentProperties& component =
            components_[i];

        const double tc =
            component.criticalTemperatureK;

        const double pc =
            component.criticalPressurePa;

        if (tc <= 0.0 || pc <= 0.0) {
            throw std::runtime_error(
                "DynamicPengRobinsonEOS: invalid critical properties"
                );
        }

        const double a0 =
            OmegaA
            * GasConstantJPerKmolK
            * GasConstantJPerKmolK
            * tc
            * tc
            / pc;

        const double tr =
            temperatureK / tc;

        const double sqrtTr =
            std::sqrt(tr);

        const double alphaRoot =
            1.0 + mPure_[i] * (1.0 - sqrtTr);

        const double alpha =
            alphaRoot * alphaRoot;

        aPure[i] =
            a0 * alpha;

        const double dSqrtTrDt =
            1.0 / (2.0 * std::sqrt(temperatureK * tc));

        const double dAlphaRootDt =
            -mPure_[i] * dSqrtTrDt;

        const double dAlphaDt =
            2.0 * alphaRoot * dAlphaRootDt;

        daPureDt[i] =
            a0 * dAlphaDt;
    }

    double aMix = 0.0;
    double daMixDt = 0.0;
    double bMix = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        bMix +=
            x[i] * bPure_[i];

        for (std::size_t j = 0; j < n; ++j) {
            const double aij =
                std::sqrt(aPure[i] * aPure[j])
                * (1.0 - kij_[i][j]);

            double daijDt = 0.0;

            if (aPure[i] > 0.0 && aPure[j] > 0.0) {
                daijDt =
                    0.5
                    * aij
                    * (
                        daPureDt[i] / aPure[i]
                        + daPureDt[j] / aPure[j]
                        );
            }

            const double coefficient =
                x[i] * x[j];

            aMix +=
                coefficient * aij;

            daMixDt +=
                coefficient * daijDt;
        }
    }

    if (bMix <= 0.0) {
        throw std::runtime_error(
            "DynamicPengRobinsonEOS: invalid bMix in departure enthalpy"
            );
    }

    const double B =
        bMix * pressurePa
        / (GasConstantJPerKmolK * temperatureK);

    const double numerator =
        zFactor + (1.0 + Sqrt2) * B;

    const double denominator =
        zFactor + (1.0 - Sqrt2) * B;

    if (numerator <= 0.0 || denominator <= 0.0) {
        throw std::runtime_error(
            "DynamicPengRobinsonEOS: invalid logarithm argument in departure enthalpy"
            );
    }

    const double logarithm =
        std::log(numerator / denominator);

    const double departure =
        GasConstantJPerKmolK
            * temperatureK
            * (zFactor - 1.0)
        + (
              temperatureK * daMixDt
              - aMix
              )
              / (2.0 * Sqrt2 * bMix)
              * logarithm;

    return departure;
}

const ComponentProperties& DynamicPengRobinsonEOS::component(
    std::size_t componentIndex
    ) const {
    return components_[componentIndex];
}

double DynamicPengRobinsonEOS::computeZFactorFromMolarVolume(
    double pressurePa,
    double temperatureK,
    double molarVolumeM3PerKmol
    ) const {
    if (pressurePa <= 0.0 || temperatureK <= 0.0 ||
        molarVolumeM3PerKmol <= 0.0) {
        throw std::runtime_error(
            "DynamicPengRobinsonEOS: invalid input for Z from molar volume"
            );
    }

    return pressurePa
           * molarVolumeM3PerKmol
           / (GasConstantJPerKmolK * temperatureK);
}

double DynamicPengRobinsonEOS::computePressureFromMolarVolume(
    double temperatureK,
    const DynamicComposition& composition,
    double molarVolumeM3PerKmol
    ) const {
    validateComposition(composition);

    if (temperatureK <= 0.0 || molarVolumeM3PerKmol <= 0.0) {
        throw std::runtime_error(
            "DynamicPengRobinsonEOS: invalid input for pressure from molar volume"
            );
    }

    const DynamicComposition x =
        normalizedComposition(composition);

    const std::size_t n =
        components_.size();

    std::vector<double> aPure(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        const ComponentProperties& c =
            components_[i];

        const double tc =
            c.criticalTemperatureK;

        const double pc =
            c.criticalPressurePa;

        if (tc <= 0.0 || pc <= 0.0) {
            throw std::runtime_error(
                "DynamicPengRobinsonEOS: invalid critical properties"
                );
        }

        const double a0 =
            OmegaA
            * GasConstantJPerKmolK
            * GasConstantJPerKmolK
            * tc
            * tc
            / pc;

        const double tr =
            temperatureK / tc;

        const double alphaRoot =
            1.0 + mPure_[i] * (1.0 - std::sqrt(tr));

        const double alpha =
            alphaRoot * alphaRoot;

        aPure[i] =
            a0 * alpha;
    }

    double aMix = 0.0;
    double bMix = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        bMix +=
            x[i] * bPure_[i];

        for (std::size_t j = 0; j < n; ++j) {
            const double aij =
                std::sqrt(aPure[i] * aPure[j])
                * (1.0 - kij_[i][j]);

            aMix +=
                x[i] * x[j] * aij;
        }
    }

    const double v =
        molarVolumeM3PerKmol;

    const double repulsionDenominator =
        v - bMix;

    const double attractionDenominator =
        v * v + 2.0 * bMix * v - bMix * bMix;

    if (repulsionDenominator <= 0.0 ||
        attractionDenominator <= 0.0) {
        throw std::runtime_error(
            "DynamicPengRobinsonEOS: invalid molar volume for PR pressure"
            );
    }

    return
        GasConstantJPerKmolK * temperatureK / repulsionDenominator
        -
        aMix / attractionDenominator;
}
