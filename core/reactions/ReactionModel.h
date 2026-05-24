#pragma once

#include "../Component.h"
#include "../thermo/Material.h"
#include "../thermo/PengRobinsonEOS.h"

#include <array>
#include <string>

struct ReactionRates
{
    // Скорости образования компонентов, kmol/(m3_liquid*s)
    Composition componentRatesKmolPerM3S = makeComposition();

    // Для диагностики
    double forwardRateKmolPerM3S = 0.0;
    double reverseRateKmolPerM3S = 0.0;
};

class ReactionModel
{
public:
    ReactionModel(
        const MaterialList& materials
        );

    ReactionRates computeRates(
        double temperatureK,
        double pressurePa,
        const Composition& liquidComposition,
        double liquidMolarDensityKmolPerM3
        ) const;

    void setForwardRateConstant(
        double value
        );

    void setReverseRateConstant(
        double value
        );

    double forwardRateConstant() const;
    double reverseRateConstant() const;
    void setEnabled(bool enabled);
    bool enabled() const;

private:
    MaterialList materials_;

    // Пока как в Aspen: k(1), k(2)
    // Единицы: 1/s, если C задана в kmol/m3,
    // тогда r = k*C имеет kmol/(m3*s)
    double kForwardPerS_ = 100.0;
    double kReversePerS_ = 100.0;
    bool enabled_ = true;
};
