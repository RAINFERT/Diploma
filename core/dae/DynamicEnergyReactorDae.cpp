#include "DynamicEnergyReactorDae.h"
#include "../models/DynamicWellStirredReactorModel.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace {

constexpr double GasConstantJPerKmolK = 8314.46261815324;

double clampValue(
    double value,
    double minValue,
    double maxValue
    ) {
    return std::max(minValue, std::min(value, maxValue));
}

} // namespace

DynamicEnergyReactorDae::DynamicEnergyReactorDae(
    const ComponentSet& components,
    const DynamicGibbsOdeSystem& gibbsOde,
    const DynamicEnthalpyModel& enthalpyModel,
    double reactorVolumeM3
    )
    : components_(components),
    gibbsOde_(gibbsOde),
    enthalpyModel_(enthalpyModel),
    layout_(components.size()),
    reactorVolumeM3_(reactorVolumeM3)
{
    if (reactorVolumeM3_ <= 0.0) {
        throw std::runtime_error(
            "DynamicEnergyReactorDae: reactor volume must be positive"
            );
    }
}
const DynamicStateLayout& DynamicEnergyReactorDae::layout() const {
    return layout_;
}

void DynamicEnergyReactorDae::checkVectorSize(
    const DynamicVector& vector,
    const char* vectorName
    ) const {
    if (vector.size() != layout_.variableCount()) {
        throw std::runtime_error(
            std::string("DynamicEnergyReactorDae: invalid size of ") +
            vectorName
            );
    }
}

void DynamicEnergyReactorDae::checkCompositionSize(
    const DynamicComposition& composition,
    const char* compositionName
    ) const {
    if (composition.size() != layout_.componentCount()) {
        throw std::runtime_error(
            std::string("DynamicEnergyReactorDae: invalid size of ") +
            compositionName
            );
    }
}

DynamicVector DynamicEnergyReactorDae::zeroVector() const {
    return DynamicVector(layout_.variableCount(), 0.0);
}

DynamicVector DynamicEnergyReactorDae::variablesFromState(
    const DynamicComposition& massesKg,
    double energyJ,
    double pressurePa,
    double temperatureK,
    const DynamicComposition& liquidSplit,
    double liquidMolarVolumeM3PerKmol,
    double vaporMolarVolumeM3PerKmol
    ) const {
    checkCompositionSize(massesKg, "massesKg");
    checkCompositionSize(liquidSplit, "liquidSplit");

    if (liquidMolarVolumeM3PerKmol <= 0.0 ||
        vaporMolarVolumeM3PerKmol <= 0.0) {
        throw std::runtime_error(
            "DynamicEnergyReactorDae: molar volumes must be positive"
            );
    }

    DynamicVector variables(layout_.variableCount(), 0.0);

    for (std::size_t i = 0; i < layout_.componentCount(); ++i) {
        variables[layout_.mass(i)] =
            massesKg[i];

        variables[layout_.liquidSplit(i)] =
            clampValue(
                liquidSplit[i],
                1.0e-12,
                1.0 - 1.0e-12
                );
    }

    variables[layout_.energy()] = energyJ;
    variables[layout_.pressure()] = pressurePa;
    variables[layout_.temperature()] = temperatureK;

    variables[layout_.liquidMolarVolume()] =
        liquidMolarVolumeM3PerKmol;

    variables[layout_.vaporMolarVolume()] =
        vaporMolarVolumeM3PerKmol;

    return variables;
}

DynamicComposition DynamicEnergyReactorDae::extractMassesKg(
    const DynamicVector& variables
    ) const {
    checkVectorSize(variables, "variables");

    DynamicComposition massesKg(layout_.componentCount(), 0.0);

    constexpr double NegativeMassToleranceKg = 1.0e-8;

    for (std::size_t i = 0; i < layout_.componentCount(); ++i) {
        massesKg[i] = variables[layout_.mass(i)];

        if (massesKg[i] < -NegativeMassToleranceKg) {
            throw std::runtime_error(
                "DynamicEnergyReactorDae: negative component mass"
                );
        }

        if (massesKg[i] < 0.0) {
            massesKg[i] = 0.0;
        }
    }

    return massesKg;
}

DynamicComposition DynamicEnergyReactorDae::extractLiquidSplit(
    const DynamicVector& variables
    ) const {
    checkVectorSize(variables, "variables");

    DynamicComposition split(layout_.componentCount(), 0.0);

    for (std::size_t i = 0; i < layout_.componentCount(); ++i) {
        split[i] =
            clampValue(
                variables[layout_.liquidSplit(i)],
                1.0e-12,
                1.0 - 1.0e-12
                );
    }

    return split;
}

double DynamicEnergyReactorDae::extractLiquidMolarVolume(
    const DynamicVector& variables
    ) const {
    checkVectorSize(variables, "variables");

    const double value =
        variables[layout_.liquidMolarVolume()];

    if (value <= 0.0) {
        throw std::runtime_error(
            "DynamicEnergyReactorDae: liquid molar volume must be positive"
            );
    }

    return value;
}

double DynamicEnergyReactorDae::extractVaporMolarVolume(
    const DynamicVector& variables
    ) const {
    checkVectorSize(variables, "variables");

    const double value =
        variables[layout_.vaporMolarVolume()];

    if (value <= 0.0) {
        throw std::runtime_error(
            "DynamicEnergyReactorDae: vapor molar volume must be positive"
            );
    }

    return value;
}

DynamicComposition DynamicEnergyReactorDae::massesToMolesKmol(
    const DynamicComposition& massesKg
    ) const {
    checkCompositionSize(massesKg, "massesKg");

    DynamicComposition molesKmol(layout_.componentCount(), 0.0);

    for (std::size_t i = 0; i < layout_.componentCount(); ++i) {
        const double molecularWeight =
            components_[i].molarMassKgPerKmol;

        if (molecularWeight <= 0.0) {
            throw std::runtime_error(
                "DynamicEnergyReactorDae: invalid molar mass"
                );
        }

        molesKmol[i] = massesKg[i] / molecularWeight;
    }

    return molesKmol;
}

DynamicComposition DynamicEnergyReactorDae::molesToOverallComposition(
    const DynamicComposition& molesKmol
    ) const {
    checkCompositionSize(molesKmol, "molesKmol");

    const double totalMoles =
        std::accumulate(
            molesKmol.begin(),
            molesKmol.end(),
            0.0
            );

    if (totalMoles <= 0.0) {
        throw std::runtime_error(
            "DynamicEnergyReactorDae: total amount must be positive"
            );
    }

    DynamicComposition z(layout_.componentCount(), 0.0);

    for (std::size_t i = 0; i < layout_.componentCount(); ++i) {
        z[i] = molesKmol[i] / totalMoles;
    }

    return z;
}

DynamicEnergyReactorDaeEvaluation DynamicEnergyReactorDae::evaluate(
    const DynamicVector& variables
    ) const {
    checkVectorSize(variables, "variables");

    DynamicEnergyReactorDaeEvaluation evaluation;

    evaluation.massesKg =
        extractMassesKg(variables);

    evaluation.molesKmol =
        massesToMolesKmol(evaluation.massesKg);

    evaluation.totalMolesKmol =
        std::accumulate(
            evaluation.molesKmol.begin(),
            evaluation.molesKmol.end(),
            0.0
            );

    evaluation.overallComposition =
        molesToOverallComposition(evaluation.molesKmol);

    evaluation.energyJ =
        variables[layout_.energy()];

    evaluation.pressurePa =
        variables[layout_.pressure()];

    evaluation.temperatureK =
        variables[layout_.temperature()];

    if (evaluation.pressurePa <= 0.0 ||
        evaluation.temperatureK <= 0.0) {
        throw std::runtime_error(
            "DynamicEnergyReactorDae: pressure and temperature must be positive"
            );
    }

    evaluation.liquidSplit =
        extractLiquidSplit(variables);

    evaluation.liquidSplit =
        extractLiquidSplit(variables);

    evaluation.gibbsOde =
        gibbsOde_.evaluate(
            evaluation.pressurePa,
            evaluation.temperatureK,
            evaluation.totalMolesKmol,
            evaluation.overallComposition,
            evaluation.liquidSplit,
            extractLiquidMolarVolume(variables),
            extractVaporMolarVolume(variables)
            );

    evaluation.calculatedVolumeM3 =
        evaluation.totalMolesKmol
        * (
            evaluation.gibbsOde.betaLiquid
                * evaluation.gibbsOde.liquidMolarVolumeM3PerKmol
            +
            evaluation.gibbsOde.betaVapor
                * evaluation.gibbsOde.vaporMolarVolumeM3PerKmol
            );

    evaluation.volumeResidual =
        (
            evaluation.calculatedVolumeM3
            - reactorVolumeM3_
            ) / reactorVolumeM3_;

    evaluation.liquidIdealEnthalpyJPerKmol =
        enthalpyModel_.vaporMixtureMolarEnthalpyJPerKmol(
            evaluation.temperatureK,
            evaluation.gibbsOde.xLiquid
            );

    evaluation.vaporIdealEnthalpyJPerKmol =
        enthalpyModel_.vaporMixtureMolarEnthalpyJPerKmol(
            evaluation.temperatureK,
            evaluation.gibbsOde.yVapor
            );

    evaluation.liquidDepartureEnthalpyJPerKmol =
        enthalpyModel_.phaseDepartureEnthalpyJPerKmol(
            evaluation.pressurePa,
            evaluation.temperatureK,
            evaluation.gibbsOde.xLiquid,
            evaluation.gibbsOde.liquidZ
            );

    evaluation.vaporDepartureEnthalpyJPerKmol =
        enthalpyModel_.phaseDepartureEnthalpyJPerKmol(
            evaluation.pressurePa,
            evaluation.temperatureK,
            evaluation.gibbsOde.yVapor,
            evaluation.gibbsOde.vaporZ
            );

    evaluation.liquidTotalEnthalpyJPerKmol =
        evaluation.liquidIdealEnthalpyJPerKmol
        + evaluation.liquidDepartureEnthalpyJPerKmol;

    evaluation.vaporTotalEnthalpyJPerKmol =
        evaluation.vaporIdealEnthalpyJPerKmol
        + evaluation.vaporDepartureEnthalpyJPerKmol;

    evaluation.twoPhaseMolarEnthalpyJPerKmol =
        evaluation.gibbsOde.betaLiquid
            * evaluation.liquidTotalEnthalpyJPerKmol
        +
        evaluation.gibbsOde.betaVapor
            * evaluation.vaporTotalEnthalpyJPerKmol;

    evaluation.calculatedInventoryEnergyJ =
        evaluation.totalMolesKmol
        * evaluation.twoPhaseMolarEnthalpyJPerKmol;

    const double energyScale =
        std::max(
            1.0,
            std::abs(evaluation.calculatedInventoryEnergyJ)
                + evaluation.totalMolesKmol * 1.0e8
            );

    evaluation.inventoryEnergyResidual =
        (
            evaluation.energyJ
            - evaluation.calculatedInventoryEnergyJ
            ) / energyScale;

    return evaluation;
}

DynamicVector DynamicEnergyReactorDae::computeResidualVector(
    const DynamicVector& variables,
    const DynamicVector& derivatives,
    const DynamicComposition& massRatesKgPerS,
    double energyRateW
    ) const {
    checkVectorSize(variables, "variables");
    checkVectorSize(derivatives, "derivatives");
    checkCompositionSize(massRatesKgPerS, "massRatesKgPerS");

    const DynamicEnergyReactorDaeEvaluation evaluation =
        evaluate(variables);

    DynamicVector residual(layout_.residualCount(), 0.0);

    for (std::size_t i = 0; i < layout_.componentCount(); ++i) {
        residual[layout_.massResidual(i)] =
            derivatives[layout_.mass(i)]
            - massRatesKgPerS[i];
    }

    residual[layout_.energyResidual()] =
        derivatives[layout_.energy()]
        - energyRateW;

    residual[layout_.volumeResidual()] =
        evaluation.volumeResidual;

    residual[layout_.inventoryEnergyResidual()] =
        evaluation.inventoryEnergyResidual;

    constexpr double MinActiveMoleFraction = 1.0e-12;
    constexpr double SinglePhaseTolerance = 1.0e-5;
    constexpr double BoundarySplitEps = 1.0e-8;
    constexpr double GasConstantJPerKmolK = 8314.46261815324;

    const bool singleLiquid =
        evaluation.gibbsOde.betaVapor <= SinglePhaseTolerance;

    const bool singleVapor =
        evaluation.gibbsOde.betaLiquid <= SinglePhaseTolerance;

    if (singleLiquid) {
        // Однофазная жидкость:
        // паровой фазы физически нет, поэтому fL=fV не требуем.
        // Фиксируем N_i^L ≈ N_i^0.
        for (std::size_t i = 0; i < layout_.componentCount(); ++i) {
            const std::size_t splitIndex =
                layout_.liquidSplit(i);

            const std::size_t residualIndex =
                layout_.liquidSplitOdeResidual(i);

            if (evaluation.overallComposition[i] <= MinActiveMoleFraction) {
                residual[residualIndex] =
                    variables[splitIndex] - BoundarySplitEps;
            } else {
                residual[residualIndex] =
                    variables[splitIndex] - (1.0 - BoundarySplitEps);
            }
        }

        // Для существующей жидкой фазы решаем EOS-условие:
        // P_EOS(T, vL, x) = P.
        residual[layout_.liquidMolarVolumeOdeResidual()] =
            (
                evaluation.gibbsOde.liquidPressureFromEosPa
                - evaluation.pressurePa
                ) / evaluation.pressurePa;

        // Паровой объем в однофазной жидкости не имеет физического смысла.
        // Фиксируем его к идеальногазовому служебному значению,
        // чтобы Newton не получил свободную переменную.
        const double idealVaporVolume =
            GasConstantJPerKmolK
            * evaluation.temperatureK
            / evaluation.pressurePa;

        residual[layout_.vaporMolarVolumeOdeResidual()] =
            (
                variables[layout_.vaporMolarVolume()]
                - idealVaporVolume
                ) / idealVaporVolume;

    } else if (singleVapor) {
        // Однофазный пар:
        // жидкой фазы физически нет, поэтому fL=fV не требуем.
        // Фиксируем N_i^L ≈ 0.
        for (std::size_t i = 0; i < layout_.componentCount(); ++i) {
            const std::size_t splitIndex =
                layout_.liquidSplit(i);

            const std::size_t residualIndex =
                layout_.liquidSplitOdeResidual(i);

            residual[residualIndex] =
                variables[splitIndex] - BoundarySplitEps;
        }

        // Жидкий объем в однофазном паре служебный.
        // Фиксируем его к текущему значению через простую положительную шкалу.
        const double dummyLiquidVolume = 1.0e-3;

        residual[layout_.liquidMolarVolumeOdeResidual()] =
            (
                variables[layout_.liquidMolarVolume()]
                - dummyLiquidVolume
                ) / dummyLiquidVolume;

        // Для существующей паровой фазы решаем EOS-условие:
        // P_EOS(T, vV, y) = P.
        residual[layout_.vaporMolarVolumeOdeResidual()] =
            (
                evaluation.gibbsOde.vaporPressureFromEosPa
                - evaluation.pressurePa
                ) / evaluation.pressurePa;

    } else {
        // Двухфазная область:
        // здесь действительно должно выполняться f_i^L = f_i^V.
        for (std::size_t i = 0; i < layout_.componentCount(); ++i) {
            const std::size_t residualIndex =
                layout_.liquidSplitOdeResidual(i);

            if (evaluation.overallComposition[i] <= MinActiveMoleFraction) {
                residual[residualIndex] =
                    variables[layout_.liquidSplit(i)] - BoundarySplitEps;
            } else {
                residual[residualIndex] =
                    evaluation.gibbsOde.fugacityLogResiduals[i];
            }
        }

        residual[layout_.liquidMolarVolumeOdeResidual()] =
            (
                evaluation.gibbsOde.liquidPressureFromEosPa
                - evaluation.pressurePa
                ) / evaluation.pressurePa;

        residual[layout_.vaporMolarVolumeOdeResidual()] =
            (
                evaluation.gibbsOde.vaporPressureFromEosPa
                - evaluation.pressurePa
                ) / evaluation.pressurePa;
    }

    return residual;
}

DynamicVector DynamicEnergyReactorDae::computeCoupledResidualVector(
    const DynamicVector& variables,
    const DynamicVector& derivatives,
    const DynamicWellStirredReactorModel& reactorModel
    ) const {
    checkVectorSize(variables, "variables");
    checkVectorSize(derivatives, "derivatives");

    const DynamicEnergyReactorDaeEvaluation evaluation =
        evaluate(variables);

    const DynamicMassBalanceResult massBalance =
        reactorModel.computeMassBalance(
            evaluation.pressurePa,
            evaluation.massesKg,
            evaluation.gibbsOde.yVapor
            );

    const DynamicEnergyBalanceResult energyBalance =
        reactorModel.computeEnergyBalance(
            evaluation.pressurePa,
            evaluation.temperatureK,
            massBalance.outletMolarFlowKmolPerS,
            evaluation.gibbsOde.yVapor
            );

    return computeResidualVector(
        variables,
        derivatives,
        massBalance.massRatesKgPerS,
        energyBalance.energyRateW
        );
}
