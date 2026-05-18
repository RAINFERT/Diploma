#include "Material.h"

#include <cmath>

MaterialList createDefaultMaterials()
{
    return {
        Material{
            Component::C2H6,
            "C2H6",
            32.27800903 + 273.15,     // Tc, K
            4.883850098e6,            // Pc, Pa
            9.86e-2,                  // omega
            30.06990051,              // kg/kmol
            AntoineCoefficients{
                44.0103,
                -2568.82,
                0.0,
                -4.97635,
                1.46e-5,
                2.0
            }
        },
        Material{
            Component::C5H12,
            "C5H12",
            196.4500061 + 273.15,     // Tc, K
            3.375120117e6,            // Pc, Pa
            0.253890008,              // omega
            72.15100098,              // kg/kmol
            AntoineCoefficients{
                63.3315,
                -5117.78,
                0.0,
                -7.48305,
                7.77e-6,
                2.0
            }
        },
        Material{
            Component::H2O,
            "H2O",
            374.1490112 + 273.15,     // Tc, K
            22.12e6,                  // Pc, Pa
            0.344000012,              // omega
            18.01510048,              // kg/kmol
            AntoineCoefficients{
                65.9278,
                -7227.53,
                0.0,
                -7.17695,
                4.03e-6,
                2.0
            }
        }
    };
}

double vaporPressurePa(const Material& material, double temperatureK)
{
    const AntoineCoefficients& c = material.antoine;

    const double lnPressureKPa =
        c.A
        + c.B / (c.C + temperatureK)
        + c.D * std::log(temperatureK)
        + c.E * std::pow(temperatureK, c.F);

    const double pressureKPa = std::exp(lnPressureKPa);

    return pressureKPa * 1000.0;
}
