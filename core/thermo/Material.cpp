#include "Material.h"

#include <cmath>
#include <stdexcept>

MaterialList createDefaultMaterials()
{
    return {
        Material{
            Component::C2H6,
            "C2H6",

            // Данные из chemsep_components_model_ready.csv:
            // component_name = Ethane
            305.32,       // Tc, K
            4872000.0,    // Pc, Pa
            0.1455,       // Vc, m3/kmol
            0.279,        // Zc
            0.099,        // omega
            30.06904,     // MW, kg/kmol

            // ChemSep vapour_pressure_eqno = 101, unit = Pa
            // ln(P[Pa]) = A + B/T + C*ln(T) + D*T^E
            AntoineCoefficients{
                61.43744,     // A
                -2814.319,    // B
                0.0,          // denominator shift, not used for eq101
                -6.778053,    // C in ChemSep, coefficient before ln(T)
                2.10827e-05,  // D
                2.0           // E
            }
        },

        Material{
            Component::C5H12,
            "C5H12",

            // Данные из chemsep_components_model_ready.csv:
            // component_name = N-pentane
            469.7,        // Tc, K
            3370000.0,    // Pc, Pa
            0.311,        // Vc, m3/kmol
            0.268,        // Zc
            0.251,        // omega
            72.14878,     // MW, kg/kmol

            // ChemSep vapour_pressure_eqno = 101, unit = Pa
            AntoineCoefficients{
                72.14242,       // A
                -5265.589,      // B
                0.0,            // denominator shift
                -7.720709,      // C in ChemSep
                7.151866e-06,   // D
                2.0             // E
            }
        },

        Material{
            Component::H2O,
            "H2O",

            // Данные из chemsep_components_model_ready.csv:
            // component_name = Water
            647.14,       // Tc, K
            22064000.0,   // Pc, Pa
            0.05595,      // Vc, m3/kmol
            0.229,        // Zc
            0.344,        // omega
            18.01528,     // MW, kg/kmol

            // ChemSep vapour_pressure_eqno = 101, unit = Pa
            AntoineCoefficients{
                74.55502,     // A
                -7295.586,    // B
                0.0,          // denominator shift
                -7.442448,    // C in ChemSep
                4.2881e-06,   // D
                2.0           // E
            }
        }
    };
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
