#include "PengRobinsonEOS.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{
constexpr double Pi = 3.14159265358979323846;
constexpr double RootTolerance = 1.0e-10;

std::vector<double> realRootsOfCubic(
    double a,
    double b,
    double c
    )
{
    // Решаем кубическое уравнение:
    //
    // x^3 + a*x^2 + b*x + c = 0
    //
    // Сначала делаем замену:
    //
    // x = u - a/3
    //
    // После этого получаем приведенное кубическое уравнение:
    //
    // u^3 + p*u + q = 0

    const double p = b - a * a / 3.0;
    const double q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;

    const double discriminant =
        q * q / 4.0 + p * p * p / 27.0;

    std::vector<double> roots;

    if (discriminant > RootTolerance)
    {
        // Один действительный корень
        const double sqrtD = std::sqrt(discriminant);

        const double u =
            std::cbrt(-q / 2.0 + sqrtD)
            + std::cbrt(-q / 2.0 - sqrtD);

        roots.push_back(u - a / 3.0);
    }
    else if (std::abs(discriminant) <= RootTolerance)
    {
        // Кратные действительные корни
        const double u1 = 2.0 * std::cbrt(-q / 2.0);
        const double u2 = -std::cbrt(-q / 2.0);

        roots.push_back(u1 - a / 3.0);
        roots.push_back(u2 - a / 3.0);
    }
    else
    {
        // Три действительных корня
        const double radius = 2.0 * std::sqrt(-p / 3.0);
        const double argument =
            (3.0 * q / (2.0 * p)) * std::sqrt(-3.0 / p);

        const double safeArgument =
            std::max(-1.0, std::min(1.0, argument));

        const double phi = std::acos(safeArgument);

        for (int k = 0; k < 3; ++k)
        {
            const double u =
                radius * std::cos((phi - 2.0 * Pi * k) / 3.0);

            roots.push_back(u - a / 3.0);
        }
    }

    std::sort(roots.begin(), roots.end());

    return roots;
}

std::vector<double> positiveUniqueRoots(const std::vector<double>& roots)
{
    std::vector<double> result;

    for (double root : roots)
    {
        if (root <= 0.0)
        {
            continue;
        }

        if (result.empty() || std::abs(root - result.back()) > RootTolerance)
        {
            result.push_back(root);
        }
    }

    return result;
}
}

PengRobinsonEOS::PengRobinsonEOS(const MaterialList& materials)
    : materials_(materials),
    bi_{},
    mi_{},
    kij_{}
{
    computePureParameters();
}

void PengRobinsonEOS::computePureParameters()
{
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const Material& material = materials_[i];

        bi_[i] =
            0.07780
            * R
            * material.criticalTemperatureK
            / material.criticalPressurePa;

        const double omega = material.acentricFactor;

        if (omega <= 0.491)
        {
            mi_[i] =
                0.37464
                + 1.54226 * omega
                - 0.26992 * omega * omega;
        }
        else
        {
            mi_[i] =
                0.379642
                + 1.48503 * omega
                - 0.164423 * omega * omega
                + 0.016666 * omega * omega * omega;
        }
    }
}

double PengRobinsonEOS::computeAlpha(
    std::size_t componentIndex,
    double temperatureK
    ) const
{
    const Material& material = materials_[componentIndex];

    const double reducedTemperature =
        temperatureK / material.criticalTemperatureK;

    const double sqrtTr = std::sqrt(reducedTemperature);

    const double value =
        1.0 + mi_[componentIndex] * (1.0 - sqrtTr);

    return value * value;
}

std::array<double, ComponentCount> PengRobinsonEOS::computeAi(
    double temperatureK
    ) const
{
    std::array<double, ComponentCount> ai{};

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const Material& material = materials_[i];

        const double alpha = computeAlpha(i, temperatureK);

        ai[i] =
            0.45724
            * R
            * R
            * material.criticalTemperatureK
            * material.criticalTemperatureK
            / material.criticalPressurePa
            * alpha;
    }

    return ai;
}

MixtureParameters PengRobinsonEOS::computeMixtureParameters(
    const Composition& composition,
    double temperatureK
    ) const
{
    const std::array<double, ComponentCount> ai = computeAi(temperatureK);

    double aMix = 0.0;
    double bMix = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        bMix += composition[i] * bi_[i];
    }

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        for (std::size_t j = 0; j < ComponentCount; ++j)
        {
            aMix +=
                composition[i]
                * composition[j]
                * std::sqrt(ai[i] * ai[j])
                * (1.0 - kij_[i][j]);
        }
    }

    return MixtureParameters{aMix, bMix};
}

ZFactors PengRobinsonEOS::computeZFactors(
    double pressurePa,
    double temperatureK,
    const Composition& composition
    ) const
{
    if (pressurePa <= 0.0)
    {
        throw std::invalid_argument("Pressure must be positive");
    }

    if (temperatureK <= 0.0)
    {
        throw std::invalid_argument("Temperature must be positive");
    }

    const MixtureParameters params =
        computeMixtureParameters(composition, temperatureK);

    const double A =
        params.a
        * pressurePa
        / (R * R * temperatureK * temperatureK);

    const double B =
        params.b
        * pressurePa
        / (R * temperatureK);

    // Кубическое уравнение Peng-Robinson:
    //
    // Z^3 - (1 - B) Z^2
    // + (A - 2B - 3B^2) Z
    // - (AB - B^2 - B^3) = 0
    //
    // Для функции realRootsOfCubic нужно привести его к виду:
    //
    // Z^3 + c2 Z^2 + c1 Z + c0 = 0

    const double c2 = -(1.0 - B);
    const double c1 = A - 2.0 * B - 3.0 * B * B;
    const double c0 = -(A * B - B * B - B * B * B);

    const std::vector<double> roots =
        positiveUniqueRoots(realRootsOfCubic(c2, c1, c0));

    if (roots.empty())
    {
        throw std::runtime_error("Peng-Robinson EOS produced no positive real roots");
    }

    if (roots.size() == 1)
    {
        return ZFactors{
            roots[0],
            roots[0],
            roots[0],
            1
        };
    }

    if (roots.size() == 2)
    {
        return ZFactors{
            roots[0],
            roots[0],
            roots[1],
            2
        };
    }

    return ZFactors{
        roots[0],
        roots[1],
        roots[2],
        3
    };
}

Composition PengRobinsonEOS::computeFugacityCoefficients(
    double pressurePa,
    double temperatureK,
    const Composition& composition,
    double zFactor
    ) const
{
    if (pressurePa <= 0.0)
    {
        throw std::invalid_argument("Pressure must be positive");
    }

    if (temperatureK <= 0.0)
    {
        throw std::invalid_argument("Temperature must be positive");
    }

    const MixtureParameters params =
        computeMixtureParameters(composition, temperatureK);

    const std::array<double, ComponentCount> ai =
        computeAi(temperatureK);

    const double A =
        params.a
        * pressurePa
        / (R * R * temperatureK * temperatureK);

    const double B =
        params.b
        * pressurePa
        / (R * temperatureK);

    if (params.a <= 0.0)
    {
        throw std::runtime_error("Invalid Peng-Robinson mixture parameter a");
    }

    if (params.b <= 0.0)
    {
        throw std::runtime_error("Invalid Peng-Robinson mixture parameter b");
    }

    if (zFactor <= B)
    {
        throw std::runtime_error("Invalid Z factor: Z must be greater than B");
    }

    Composition phi{};

    const double sqrt2 = std::sqrt(2.0);

    const double logTerm =
        std::log(
            (zFactor + (1.0 + sqrt2) * B)
            / (zFactor + (1.0 - sqrt2) * B)
            );

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        double daDni = 0.0;

        for (std::size_t j = 0; j < ComponentCount; ++j)
        {
            daDni +=
                composition[j]
                * std::sqrt(ai[i] * ai[j])
                * (1.0 - kij_[i][j]);
        }

        daDni *= 2.0;

        const double Bi =
            bi_[i]
            * pressurePa
            / (R * temperatureK);

        const double lnPhi =
            Bi / B * (zFactor - 1.0)
            - std::log(zFactor - B)
            - A / (2.0 * sqrt2 * B)
                  * (daDni / params.a - Bi / B)
                  * logTerm;

        phi[i] = std::exp(lnPhi);
    }

    return phi;
}
