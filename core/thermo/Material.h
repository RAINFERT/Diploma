#pragma once

#include "../Component.h"

#include <array>
#include <string>

struct AntoineCoefficients
{
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

    // Ацентрический фактор
    double acentricFactor;

    // Молярная масса, kg/kmol
    double molarMassKgPerKmol;

    // Коэффициенты Антуана / DIPPR-подобного уравнения
    AntoineCoefficients antoine;
};

using MaterialList = std::array<Material, ComponentCount>;

MaterialList createDefaultMaterials();

double vaporPressurePa(const Material& material, double temperatureK);
