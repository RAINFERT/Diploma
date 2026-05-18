#include "EnthalpyModel.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>

double EnthalpyModel::idealGasMixtureMolarEnthalpyJPerKmol(
    double temperatureK,
    const Composition& composition
    ) const
{
    return mixtureMolarEnthalpyJPerKmol(
        temperatureK,
        composition
        );
}

namespace
{
double binaryInteractionParameter(
    std::size_t i,
    std::size_t j
    )
{
    (void)i;
    (void)j;

    return 0.0;
}

double pengRobinsonKappa(double omega)
{
    if (omega <= 0.491)
    {
        return
            0.37464
            + 1.54226 * omega
            - 0.26992 * omega * omega;
    }

    return
        0.379642
        + 1.48503 * omega
        - 0.164423 * omega * omega
        + 0.016666 * omega * omega * omega;
}
}

double EnthalpyModel::pengRobinsonDepartureEnthalpyJPerKmol(
    double pressurePa,
    double temperatureK,
    const Composition& composition,
    double compressibilityFactor
    ) const
{
    validateTemperature(temperatureK);

    if (pressurePa <= 0.0)
    {
        throw std::invalid_argument("Pressure must be positive for PR enthalpy departure");
    }

    if (compressibilityFactor <= 0.0)
    {
        throw std::invalid_argument("Compressibility factor must be positive for PR enthalpy departure");
    }

    constexpr double R = 8314.46261815324; // J/(kmol*K)
    constexpr double sqrt2 = 1.4142135623730951;

    std::array<double, ComponentCount> aPure{};
    std::array<double, ComponentCount> daPureDt{};
    std::array<double, ComponentCount> bPure{};

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const Material& material =
            materials_[i];

        const double tc =
            material.criticalTemperatureK;

        const double pc =
            material.criticalPressurePa;

        const double omega =
            material.acentricFactor;

        if (tc <= 0.0 || pc <= 0.0)
        {
            throw std::runtime_error("Invalid critical property in PR enthalpy departure");
        }

        const double tr =
            temperatureK / tc;

        const double sqrtTr =
            std::sqrt(tr);

        const double kappa =
            pengRobinsonKappa(omega);

        const double alphaRoot =
            1.0 + kappa * (1.0 - sqrtTr);

        const double alpha =
            alphaRoot * alphaRoot;

        const double a0 =
            0.45724 * R * R * tc * tc / pc;

        aPure[i] =
            a0 * alpha;

        const double dSqrtTrDt =
            1.0 / (2.0 * std::sqrt(temperatureK * tc));

        const double dAlphaRootDt =
            -kappa * dSqrtTrDt;

        const double dAlphaDt =
            2.0 * alphaRoot * dAlphaRootDt;

        daPureDt[i] =
            a0 * dAlphaDt;

        bPure[i] =
            0.07780 * R * tc / pc;
    }

    double compositionSum = 0.0;

    for (double value : composition)
    {
        if (value < 0.0)
        {
            throw std::invalid_argument("Composition contains negative value in PR enthalpy departure");
        }

        compositionSum += value;
    }

    if (compositionSum <= 0.0)
    {
        throw std::invalid_argument("Composition sum must be positive in PR enthalpy departure");
    }

    Composition normalizedComposition{};

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        normalizedComposition[i] =
            composition[i] / compositionSum;
    }

    double aMix = 0.0;
    double daMixDt = 0.0;
    double bMix = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        bMix +=
            normalizedComposition[i] * bPure[i];

        for (std::size_t j = 0; j < ComponentCount; ++j)
        {
            const double oneMinusKij =
                1.0 - binaryInteractionParameter(i, j);

            const double aij =
                std::sqrt(aPure[i] * aPure[j])
                * oneMinusKij;

            double daijDt = 0.0;

            if (aPure[i] > 0.0 && aPure[j] > 0.0)
            {
                daijDt =
                    0.5
                    * aij
                    * (
                        daPureDt[i] / aPure[i]
                        + daPureDt[j] / aPure[j]
                        );
            }

            const double coefficient =
                normalizedComposition[i]
                * normalizedComposition[j];

            aMix +=
                coefficient * aij;

            daMixDt +=
                coefficient * daijDt;
        }
    }

    if (bMix <= 0.0)
    {
        throw std::runtime_error("Invalid bMix in PR enthalpy departure");
    }

    const double B =
        bMix * pressurePa / (R * temperatureK);

    const double numerator =
        compressibilityFactor
        + (1.0 + sqrt2) * B;

    const double denominator =
        compressibilityFactor
        + (1.0 - sqrt2) * B;

    if (numerator <= 0.0 || denominator <= 0.0)
    {
        throw std::runtime_error("Invalid logarithm argument in PR enthalpy departure");
    }

    const double logarithm =
        std::log(numerator / denominator);

    const double departure =
        R * temperatureK * (compressibilityFactor - 1.0)
        + (
              temperatureK * daMixDt
              - aMix
              )
              / (2.0 * sqrt2 * bMix)
              * logarithm;

    return departure;
}


double EnthalpyModel::phaseMolarEnthalpyJPerKmol(
    double pressurePa,
    double temperatureK,
    const Composition& composition,
    double compressibilityFactor
    ) const
{
    const double idealPart =
        idealGasMixtureMolarEnthalpyJPerKmol(
            temperatureK,
            composition
            );

    const double departurePart =
        pengRobinsonDepartureEnthalpyJPerKmol(
            pressurePa,
            temperatureK,
            composition,
            compressibilityFactor
            );

    return idealPart + departurePart;
}

ComponentEnthalpyDataList createDefaultEnthalpyData()
{
    ComponentEnthalpyDataList data{};

    // ВАЖНО:
    // Это временные приближенные Cp для проверки энергетического баланса.
    // Сейчас используем Cp(T) = cpA.
    // То есть cpB = cpC = cpD = cpE = 0.
    // Позже эти коэффициенты должны загружаться из файла/БД.

    for (ComponentEnthalpyData& item : data)
    {
        item.referenceTemperatureK = 298.15;
        item.referenceMolarEnthalpyJPerKmol = 0.0;

        item.cpA = 0.0;
        item.cpB = 0.0;
        item.cpC = 0.0;
        item.cpD = 0.0;
        item.cpE = 0.0;
    }

    // C2H6, примерно газовая Cp около 52.5 J/(mol*K)
    // 52.5 J/(mol*K) = 52.5e3 J/(kmol*K)
    data[componentIndex(Component::C2H6)].name = "C2H6";
    data[componentIndex(Component::C2H6)].referenceMolarEnthalpyJPerKmol =
        -84738.0e3;
    data[componentIndex(Component::C2H6)].cpA =
        52.5e3;

    data[componentIndex(Component::C5H12)].name = "C5H12";
    data[componentIndex(Component::C5H12)].referenceMolarEnthalpyJPerKmol =
        -146490.0e3;
    data[componentIndex(Component::C5H12)].cpA =
        167.0e3;

    data[componentIndex(Component::H2O)].name = "H2O";
    data[componentIndex(Component::H2O)].referenceMolarEnthalpyJPerKmol =
        -241814.0e3;
    data[componentIndex(Component::H2O)].cpA =
        75.3e3;

    return data;
}

EnthalpyModel::EnthalpyModel(
    const MaterialList& materials
    )
    : materials_(materials),
    enthalpyData_(createDefaultEnthalpyData())
{
}

EnthalpyModel::EnthalpyModel(
    const MaterialList& materials,
    const ComponentEnthalpyDataList& enthalpyData
    )
    : materials_(materials),
    enthalpyData_(enthalpyData)
{
}

const ComponentEnthalpyDataList& EnthalpyModel::data() const
{
    return enthalpyData_;
}

ComponentEnthalpyData EnthalpyModel::componentData(
    std::size_t componentIndexValue
    ) const
{
    if (componentIndexValue >= ComponentCount)
    {
        throw std::out_of_range("Invalid component index in enthalpy data");
    }

    return enthalpyData_[componentIndexValue];
}

void EnthalpyModel::validateTemperature(
    double temperatureK
    ) const
{
    if (!std::isfinite(temperatureK) || temperatureK <= 0.0)
    {
        throw std::invalid_argument("Temperature must be positive for enthalpy calculation");
    }
}

double EnthalpyModel::componentMolarEnthalpyJPerKmol(
    std::size_t componentIndex,
    double temperatureK
    ) const
{
    validateTemperature(temperatureK);

    const ComponentEnthalpyData& component =
        componentData(componentIndex);

    const double T =
        temperatureK;

    const double Tref =
        component.referenceTemperatureK;

    const double deltaH =
        component.cpA * (T - Tref)
        + 0.5 * component.cpB * (T * T - Tref * Tref)
        + (1.0 / 3.0) * component.cpC * (T * T * T - Tref * Tref * Tref)
        + 0.25 * component.cpD * (
              T * T * T * T
              - Tref * Tref * Tref * Tref
              )
        - component.cpE * (
              1.0 / T
              - 1.0 / Tref
              );

    return
        component.referenceMolarEnthalpyJPerKmol
        + deltaH;
}

Composition EnthalpyModel::componentMolarEnthalpiesJPerKmol(
    double temperatureK
    ) const
{
    validateTemperature(temperatureK);

    Composition result{};

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        result[i] =
            componentMolarEnthalpyJPerKmol(
                i,
                temperatureK
                );
    }

    return result;
}

double EnthalpyModel::mixtureMolarEnthalpyJPerKmol(
    double temperatureK,
    const Composition& composition
    ) const
{
    validateTemperature(temperatureK);

    double compositionSum = 0.0;
    double enthalpySum = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        if (composition[i] < 0.0)
        {
            throw std::invalid_argument("Composition contains negative value");
        }

        const double hComponent =
            componentMolarEnthalpyJPerKmol(
                i,
                temperatureK
                );

        compositionSum +=
            composition[i];

        enthalpySum +=
            composition[i] * hComponent;
    }

    if (compositionSum <= 0.0)
    {
        throw std::invalid_argument("Composition sum must be positive");
    }

    // Делим на сумму на случай, если состав не идеально нормирован.
    return enthalpySum / compositionSum;
}

double EnthalpyModel::enthalpyFlowW(
    double flowKmolPerS,
    double temperatureK,
    const Composition& composition
    ) const
{
    const double h =
        mixtureMolarEnthalpyJPerKmol(
            temperatureK,
            composition
            );

    return flowKmolPerS * h;
}
