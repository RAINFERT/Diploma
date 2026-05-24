#include "Material.h"

#include <cmath>
#include <stdexcept>

MaterialList createDefaultMaterials()
{
    MaterialList materials(ComponentCount);

    materials[componentIndex(Component::C2H6)] =
        Material{
            Component::C2H6,
            "C2H6",

            // component_name = Ethane
            305.32,       // Tc, K
            4872000.0,    // Pc, Pa
            0.1455,       // Vc, m3/kmol
            0.279,        // Zc
            0.099,        // omega
            30.06904,     // MW, kg/kmol

            AntoineCoefficients{
                61.43744,
                -2814.319,
                0.0,
                -6.778053,
                2.10827e-05,
                2.0
            }
        };

    materials[componentIndex(Component::C5H12)] =
        Material{
            Component::C5H12,
            "C5H12",

            // component_name = N-pentane
            469.7,
            3370000.0,
            0.311,
            0.268,
            0.251,
            72.14878,

            AntoineCoefficients{
                72.14242,
                -5265.589,
                0.0,
                -7.720709,
                7.151866e-06,
                2.0
            }
        };

    materials[componentIndex(Component::H2O)] =
        Material{
            Component::H2O,
            "H2O",

            // component_name = Water
            647.14,
            22064000.0,
            0.05595,
            0.229,
            0.344,
            18.01528,

            AntoineCoefficients{
                74.55502,
                -7295.586,
                0.0,
                -7.442448,
                4.2881e-06,
                2.0
            }
        };

    return materials;
}

double vaporPressurePa(
    const Material& material,
    double temperatureK
    )
{
    if (temperatureK <= 0.0)
    {
        throw std::invalid_argument("Temperature must be positive");
    }

    const AntoineCoefficients& c = material.antoine;

    // Для текущих коэффициентов из ChemSep:
    //
    // ln(P[Pa]) = A + B / T + C * ln(T) + D * T^E
    //
    // В структуре:
    // c.C = 0.0
    // c.D = коэффициент перед ln(T)
    // c.E = коэффициент перед T^F
    // c.F = степень T
    const double lnPressurePa =
        c.A
        + c.B / (c.C + temperatureK)
        + c.D * std::log(temperatureK)
        + c.E * std::pow(temperatureK, c.F);

    return std::exp(lnPressurePa);
}
