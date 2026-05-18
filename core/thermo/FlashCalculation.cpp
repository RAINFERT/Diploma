#include "FlashCalculation.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
constexpr double Tiny = 1.0e-14;

double sumComposition(const Composition& composition)
{
    double sum = 0.0;

    for (double value : composition)
    {
        sum += value;
    }

    return sum;
}

Composition normalizedComposition(const Composition& composition)
{
    const double sum = sumComposition(composition);

    if (sum <= 0.0)
    {
        throw std::runtime_error("Composition sum must be positive");
    }

    Composition result{};

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        result[i] = composition[i] / sum;
    }

    return result;
}
}

std::string flashStatusToString(FlashStatus status)
{
    switch (status)
    {
    case FlashStatus::Converged:
        return "converged";
    case FlashStatus::SingleLiquid:
        return "single_liquid";
    case FlashStatus::SingleVapor:
        return "single_vapor";
    case FlashStatus::NotConverged:
        return "not_converged";
    }

    return "unknown";
}

FlashCalculation::FlashCalculation(
    const MaterialList& materials,
    const PengRobinsonEOS& eos
    )
    : materials_(materials),
    eos_(eos)
{
}

Composition FlashCalculation::wilsonKValues(
    double pressurePa,
    double temperatureK
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

    Composition k{};

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const Material& material = materials_[i];

        k[i] =
            material.criticalPressurePa
            / pressurePa
            * std::exp(
                5.31
                * (1.0 + material.acentricFactor)
                * (1.0 - material.criticalTemperatureK / temperatureK)
                );

        if (k[i] < Tiny)
        {
            k[i] = Tiny;
        }
    }

    return k;
}

Composition FlashCalculation::initialKValues(
    double pressurePa,
    double temperatureK,
    const Composition& zOverall
    ) const
{
    Composition k = wilsonKValues(pressurePa, temperatureK);

    try
    {
        const ZFactors zFactors =
            eos_.computeZFactors(
                pressurePa,
                temperatureK,
                zOverall
                );

        if (zFactors.rootCount < 2)
        {
            return k;
        }

        const Composition phiLiquid =
            eos_.computeFugacityCoefficients(
                pressurePa,
                temperatureK,
                zOverall,
                zFactors.liquid
                );

        const Composition phiVapor =
            eos_.computeFugacityCoefficients(
                pressurePa,
                temperatureK,
                zOverall,
                zFactors.vapor
                );

        Composition kPr{};

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            if (phiVapor[i] <= Tiny)
            {
                return k;
            }

            kPr[i] = phiLiquid[i] / phiVapor[i];

            if (!std::isfinite(kPr[i]) || kPr[i] < Tiny)
            {
                return k;
            }
        }

        return kPr;
    }
    catch (...)
    {
        return k;
    }
}

double FlashCalculation::rachfordRice(
    double beta,
    const Composition& z,
    const Composition& k
    ) const
{
    double value = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const double denominator =
            1.0 + beta * (k[i] - 1.0);

        if (std::abs(denominator) < Tiny)
        {
            throw std::runtime_error("Invalid denominator in Rachford-Rice equation");
        }

        value += z[i] * (k[i] - 1.0) / denominator;
    }

    return value;
}

double FlashCalculation::solveBetaBisection(
    const Composition& z,
    const Composition& k
    ) const
{
    const double f0 = rachfordRice(0.0, z, k);
    const double f1 = rachfordRice(1.0, z, k);

    // Если F(0) <= 0, корня на интервале [0, 1] нет,
    // смесь находится в жидкой области.
    if (f0 <= 0.0)
    {
        return 0.0;
    }

    // Если F(1) >= 0, корня на интервале [0, 1] нет,
    // смесь находится в паровой области.
    if (f1 >= 0.0)
    {
        return 1.0;
    }

    double left = 0.0;
    double right = 1.0;

    for (int iteration = 0; iteration < 100; ++iteration)
    {
        const double middle = 0.5 * (left + right);
        const double fMiddle = rachfordRice(middle, z, k);

        if (std::abs(fMiddle) < 1.0e-12)
        {
            return middle;
        }

        if (fMiddle > 0.0)
        {
            left = middle;
        }
        else
        {
            right = middle;
        }
    }

    return 0.5 * (left + right);
}

double FlashCalculation::solveBeta(
    const Composition& z,
    const Composition& k
    ) const
{
    const double f0 = rachfordRice(0.0, z, k);
    const double f1 = rachfordRice(1.0, z, k);

    // Однофазная жидкая область
    if (f0 <= 0.0)
    {
        return 0.0;
    }

    // Однофазная паровая область
    if (f1 >= 0.0)
    {
        return 1.0;
    }

    double left = 0.0;
    double right = 1.0;

    // Начальное приближение.
    // Для надежности берем середину физического интервала.
    double beta = 0.5;

    constexpr int MaxNewtonIterations = 20;
    constexpr double FunctionTolerance = 1.0e-12;
    constexpr double StepTolerance = 1.0e-12;

    try
    {
        for (int iteration = 0; iteration < MaxNewtonIterations; ++iteration)
        {
            const double fBeta =
                rachfordRice(beta, z, k);

            if (std::abs(fBeta) < FunctionTolerance)
            {
                return beta;
            }

            // Обновляем интервал, в котором зажат корень.
            // Rachford-Rice F(beta) монотонно убывает.
            if (fBeta > 0.0)
            {
                left = beta;
            }
            else
            {
                right = beta;
            }

            const double derivative =
                rachfordRiceDerivative(beta, z, k);

            double betaCandidate =
                beta;

            bool useNewtonStep =
                std::isfinite(derivative)
                && std::abs(derivative) > Tiny;

            if (useNewtonStep)
            {
                betaCandidate =
                    beta - fBeta / derivative;

                // Если Newton вышел за физический интервал,
                // заменяем его на безопасный шаг бисекции.
                if (!std::isfinite(betaCandidate)
                    || betaCandidate <= left
                    || betaCandidate >= right)
                {
                    useNewtonStep = false;
                }
            }

            if (!useNewtonStep)
            {
                betaCandidate =
                    0.5 * (left + right);
            }

            if (std::abs(betaCandidate - beta) < StepTolerance)
            {
                return betaCandidate;
            }

            beta = betaCandidate;
        }
    }
    catch (...)
    {
        return solveBetaBisection(z, k);
    }

    // Если Newton не сошелся за MaxNewtonIterations,
    // возвращаемся к надежной бисекции.
    return solveBetaBisection(z, k);
}

double FlashCalculation::rachfordRiceDerivative(
    double beta,
    const Composition& z,
    const Composition& k
    ) const
{
    double derivative = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const double kMinusOne =
            k[i] - 1.0;

        const double denominator =
            1.0 + beta * kMinusOne;

        if (std::abs(denominator) < Tiny)
        {
            throw std::runtime_error("Invalid denominator in Rachford-Rice derivative");
        }

        derivative -=
            z[i]
            * kMinusOne
            * kMinusOne
            / (denominator * denominator);
    }

    return derivative;
}

void FlashCalculation::computePhaseCompositions(
    double beta,
    const Composition& z,
    const Composition& k,
    Composition& x,
    Composition& y
    ) const
{
    if (beta <= 0.0)
    {
        x = z;
        y = z;
        return;
    }

    if (beta >= 1.0)
    {
        x = z;
        y = z;
        return;
    }

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const double denominator =
            1.0 + beta * (k[i] - 1.0);

        if (denominator <= 0.0)
        {
            throw std::runtime_error("Invalid phase composition denominator");
        }

        x[i] = z[i] / denominator;
        y[i] = k[i] * x[i];
    }

    x = normalizedComposition(x);
    y = normalizedComposition(y);
}

double FlashCalculation::maxRelativeDifference(
    const Composition& a,
    const Composition& b
    ) const
{
    double maxError = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const double denominator =
            std::max(std::abs(b[i]), Tiny);

        const double error =
            std::abs(a[i] - b[i]) / denominator;

        maxError = std::max(maxError, error);
    }

    return maxError;
}

double FlashCalculation::averageMolarMassKgPerKmol(
    const Composition& composition
    ) const
{
    double molarMass = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        molarMass +=
            composition[i]
            * materials_[i].molarMassKgPerKmol;
    }

    return molarMass;
}

double FlashCalculation::molarDensityKmolPerM3(
    double pressurePa,
    double temperatureK,
    double zFactor
    ) const
{
    if (zFactor <= 0.0)
    {
        throw std::runtime_error("Z factor must be positive for density calculation");
    }

    // P / (ZRT) дает mol/m3.
    // Делим на 1000, чтобы получить kmol/m3.
    return pressurePa / (zFactor * R * temperatureK) / 1000.0;
}

double FlashCalculation::massDensityKgPerM3(
    double pressurePa,
    double temperatureK,
    const Composition& composition,
    double zFactor
    ) const
{
    const double rhoMolar =
        molarDensityKmolPerM3(
            pressurePa,
            temperatureK,
            zFactor
            );

    const double molarMass =
        averageMolarMassKgPerKmol(composition);

    // kmol/m3 * kg/kmol = kg/m3
    return rhoMolar * molarMass;
}

bool FlashCalculation::isAqueousLiquid(
    const Composition& liquidComposition
    ) const
{
    const std::size_t waterIndex =
        componentIndex(Component::H2O);

    return liquidComposition[waterIndex] > 0.95;
}

double FlashCalculation::waterMassDensityKgPerM3(
    double pressurePa,
    double temperatureK
    ) const
{
    const double temperatureC = temperatureK - 273.15;

    // Простая корреляция плотности чистой воды при атмосферном давлении.
    // Хорошо работает как первое приближение в диапазоне примерно 0...100 C.
    //
    // Давление учитываем грубо через сжимаемость воды.
    // При давлениях порядка нескольких bar поправка маленькая.

    const double numerator =
        (temperatureC + 288.9414)
        * std::pow(temperatureC - 3.9863, 2.0);

    const double denominator =
        508929.2
        * (temperatureC + 68.12963);

    double density =
        1000.0
        * (1.0 - numerator / denominator);

    const double referencePressurePa = 101325.0;
    const double waterBulkModulusPa = 2.2e9;

    density *=
        1.0
        + (pressurePa - referencePressurePa)
              / waterBulkModulusPa;

    return density;
}

void FlashCalculation::applyAqueousDensityCorrection(
    FlashResult& result,
    double pressurePa,
    double temperatureK
    ) const
{
    if (!isAqueousLiquid(result.xLiquid))
    {
        return;
    }

    result.massDensityLiquidKgPerM3 =
        waterMassDensityKgPerM3(
            pressurePa,
            temperatureK
            );

    result.molarDensityLiquidKmolPerM3 =
        result.massDensityLiquidKgPerM3
        / result.molarMassLiquidKgPerKmol;
}

void FlashCalculation::fillDensities(
    FlashResult& result,
    double pressurePa,
    double temperatureK
    ) const
{
    result.molarMassLiquidKgPerKmol =
        averageMolarMassKgPerKmol(result.xLiquid);

    result.molarMassVaporKgPerKmol =
        averageMolarMassKgPerKmol(result.yVapor);

    result.molarDensityLiquidKmolPerM3 =
        molarDensityKmolPerM3(
            pressurePa,
            temperatureK,
            result.zLiquid
            );

    result.molarDensityVaporKmolPerM3 =
        molarDensityKmolPerM3(
            pressurePa,
            temperatureK,
            result.zVapor
            );

    result.massDensityLiquidKgPerM3 =
        result.molarDensityLiquidKmolPerM3
        * result.molarMassLiquidKgPerKmol;

    result.massDensityVaporKgPerM3 =
        result.molarDensityVaporKmolPerM3
        * result.molarMassVaporKgPerKmol;

    applyAqueousDensityCorrection(
        result,
        pressurePa,
        temperatureK
        );

    if (result.status == FlashStatus::SingleLiquid)
    {
        result.molarDensityMixtureKmolPerM3 =
            result.molarDensityLiquidKmolPerM3;

        result.massDensityMixtureKgPerM3 =
            result.massDensityLiquidKgPerM3;

        return;
    }

    if (result.status == FlashStatus::SingleVapor)
    {
        result.molarDensityMixtureKmolPerM3 =
            result.molarDensityVaporKmolPerM3;

        result.massDensityMixtureKgPerM3 =
            result.massDensityVaporKgPerM3;

        return;
    }

    // beta — мольная доля пара.
    // Поэтому для молярной плотности смеси используем гармоническое смешение:
    //
    // 1 / rho_mix = beta / rho_vap + (1 - beta) / rho_liq

    result.molarDensityMixtureKmolPerM3 =
        1.0 / (
            result.beta / result.molarDensityVaporKmolPerM3
            + (1.0 - result.beta) / result.molarDensityLiquidKmolPerM3
            );

    const double mixtureMolarMass =
        (1.0 - result.beta) * result.molarMassLiquidKgPerKmol
        + result.beta * result.molarMassVaporKgPerKmol;

    result.massDensityMixtureKgPerM3 =
        result.molarDensityMixtureKmolPerM3
        * mixtureMolarMass;
}

FlashResult FlashCalculation::calculate(
    double pressurePa,
    double temperatureK,
    const Composition& zOverall,
    double tolerance,
    int maxIterations
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

    if (tolerance <= 0.0)
    {
        throw std::invalid_argument("Tolerance must be positive");
    }

    if (maxIterations <= 0)
    {
        throw std::invalid_argument("Maximum iteration count must be positive");
    }

    const Composition z = normalizedComposition(zOverall);

    FlashResult result;

    Composition k = initialKValues(
        pressurePa,
        temperatureK,
        z
        );

    for (int iteration = 0; iteration < maxIterations; ++iteration)
    {
        const double beta = solveBeta(z, k);

        Composition x{};
        Composition y{};

        computePhaseCompositions(beta, z, k, x, y);

        if (beta <= 0.0)
        {
            const ZFactors zFactors =
                eos_.computeZFactors(pressurePa, temperatureK, x);

            result.status = FlashStatus::SingleLiquid;
            result.iterations = iteration + 1;
            result.beta = 0.0;
            result.xLiquid = x;
            result.yVapor = y;
            result.kValues = k;
            result.zLiquid = zFactors.liquid;
            result.zVapor = zFactors.vapor;

            fillDensities(result, pressurePa, temperatureK);

            return result;
        }

        if (beta >= 1.0)
        {
            const ZFactors zFactors =
                eos_.computeZFactors(pressurePa, temperatureK, y);

            result.status = FlashStatus::SingleVapor;
            result.iterations = iteration + 1;
            result.beta = 1.0;
            result.xLiquid = x;
            result.yVapor = y;
            result.kValues = k;
            result.zLiquid = zFactors.liquid;
            result.zVapor = zFactors.vapor;

            fillDensities(result, pressurePa, temperatureK);

            return result;
        }

        const ZFactors zLiquidFactors =
            eos_.computeZFactors(pressurePa, temperatureK, x);

        const ZFactors zVaporFactors =
            eos_.computeZFactors(pressurePa, temperatureK, y);

        const double zLiquid = zLiquidFactors.liquid;
        const double zVapor = zVaporFactors.vapor;

        const Composition phiLiquid =
            eos_.computeFugacityCoefficients(
                pressurePa,
                temperatureK,
                x,
                zLiquid
                );

        const Composition phiVapor =
            eos_.computeFugacityCoefficients(
                pressurePa,
                temperatureK,
                y,
                zVapor
                );

        Composition kNew{};

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            kNew[i] = phiLiquid[i] / phiVapor[i];

            if (kNew[i] < Tiny)
            {
                kNew[i] = Tiny;
            }
        }

        const double error = maxRelativeDifference(kNew, k);

        result.status = FlashStatus::NotConverged;
        result.iterations = iteration + 1;
        result.beta = beta;
        result.xLiquid = x;
        result.yVapor = y;
        result.kValues = kNew;
        result.zLiquid = zLiquid;
        result.zVapor = zVapor;
        result.phiLiquid = phiLiquid;
        result.phiVapor = phiVapor;

        if (error < tolerance)
        {
            result.status = FlashStatus::Converged;
            fillDensities(result, pressurePa, temperatureK);
            return result;
        }

        // Демпфирование, как в Python-коде:
        // K = 0.1*K_old + 0.9*K_new
        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            k[i] = 0.1 * k[i] + 0.9 * kNew[i];

            if (k[i] < Tiny)
            {
                k[i] = Tiny;
            }
        }
    }

    fillDensities(result, pressurePa, temperatureK);
    return result;
}
