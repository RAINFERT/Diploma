#pragma once

#include "../Component.h"

#include <vector>
#include <string>

struct AntoineCoefficients
{
    // Используем форму:
    //
    // ln(P) = A + B / (C + T) + D * ln(T) + E * T^F
    //
    // Для ChemSep vapour_pressure_eqno = 101:
    //
    // ln(P) = A + B / T + C_chemsep * ln(T) + D_chemsep * T^E_chemsep
    //
    // Поэтому при записи из ChemSep:
    // A = vapour_pressure_A
    // B = vapour_pressure_B
    // C = 0.0
    // D = vapour_pressure_C
    // E = vapour_pressure_D
    // F = vapour_pressure_E
    double A;
    double B;
    double C;
    double D;
    double E;
    double F;
};

struct Material
{
    Component component;
    std::string name;

    // Критическая температура, K
    double criticalTemperatureK;

    // Критическое давление, Pa
    double criticalPressurePa;

    // Критический молярный объем, m3/kmol
    double criticalVolumeM3PerKmol;

    // Критический коэффициент сжимаемости
    double criticalCompressibility;

    // Ацентрический фактор
    double acentricFactor;

    // Молярная масса, kg/kmol
    double molarMassKgPerKmol;

    // Коэффициенты давления насыщения
    AntoineCoefficients antoine;
};

using MaterialList = std::vector<Material>;

MaterialList createDefaultMaterials();

double vaporPressurePa(
    const Material& material,
    double temperatureK
    );
