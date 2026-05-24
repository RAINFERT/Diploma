#pragma once

#include "../ReactorState.h"
#include "../thermo/ThermoPackage.h"
#include "../reactions/ReactionModel.h"

#include <array>
#include <string>

struct WellStirredReactorParameters
{
    // Заданный объем реактора, м3
    double volumeM3 = 1.0;

    // Входной поток, kmol/s
    double inletFlowKmolPerS = 0.001;

    // Входные условия
    double inletPressureBar = 2.5;
    double inletTemperatureC = 25.0;

    // Входной состав
    Composition inletComposition = makeLegacyComposition(0.5, 0.5, 0.0);

    // Регулятор выходного потока:
    // raw Outlet.F = outletValveCoefficient * (P - outletReferencePressure)
    double outletValveCoefficientKmolPerSBar = 0.1;
    double outletReferencePressureBar = 2.5;

    // true  — разрешить отрицательный выходной расход, как в Aspen-заглушке
    // false — запретить обратный поток: Fout = max(0, rawFout)
    bool allowReverseOutletFlow = true;

    // Пока просто храним, позже используем в энергетическом балансе
    double heatTransferAreaM2 = 1.0;

    // W/(m2*K)
    double heatTransferCoefficientWPerM2K = 500.0;

    double jacketTemperatureC = 110.0;
};

struct WellStirredReactorEvaluation
{
    // Давление и температура в SI
    double pressurePa = 0.0;
    double temperatureK = 0.0;

    // Количество вещества в реакторе
    double totalMolesKmol = 0.0;

    // Общий состав смеси
    Composition zOverall = makeComposition();

    // Результат flash-расчета
    FlashResult flash;

    // Объем, рассчитанный из количества вещества и плотности
    double calculatedVolumeM3 = 0.0;

    // Невязка уравнения объема:
    // calculatedVolumeM3 - заданный volumeM3
    double volumeResidualM3 = 0.0;
};

struct PressureInitializationResult
{
    bool converged = false;
    int iterations = 0;

    double pressurePa = 0.0;
    double pressureBar = 0.0;

    double residualM3 = 0.0;

    // Для диагностики: newton, bisection, newton_failed_then_bisection
    std::string method;

    // Сообщение от численного метода
    std::string message;

    ReactorState initializedState;
    WellStirredReactorEvaluation evaluation;
};

struct MassBalanceResult
{
    WellStirredReactorEvaluation evaluation;

    double inletFlowKmolPerS = 0.0;
    double outletFlowKmolPerS = 0.0;
    double rawOutletFlowKmolPerS = 0.0;

    Composition inletComposition = makeComposition();
    Composition outletComposition = makeComposition();

    // Объемная доля жидкой фазы внутри реактора
    double liquidVolumeFraction = 0.0;

    // Скорости реакций, kmol/(m3_liquid*s)
    ReactionRates reactionRates;

    // Отдельные вклады в dM/dt, kg/s
    Composition inletMassFlowKgPerS = makeComposition();
    Composition outletMassFlowKgPerS = makeComposition();
    Composition reactionMassRateKgPerS = makeComposition();
    Composition massDerivativesKgPerS = makeComposition();

    double totalInletMassFlowKgPerS = 0.0;
    double totalOutletMassFlowKgPerS = 0.0;
    double totalReactionMassRateKgPerS = 0.0;
    double totalMassDerivativeKgPerS = 0.0;
};

struct EnergyBalanceResult
{
    WellStirredReactorEvaluation evaluation;

    double inletMolarEnthalpyJPerKmol = 0.0;
    double outletMolarEnthalpyJPerKmol = 0.0;
    double reactorMolarEnthalpyJPerKmol = 0.0;

    double inletEnthalpyFlowW = 0.0;
    double outletEnthalpyFlowW = 0.0;

    double heatTransferW = 0.0;

    // Итоговая производная энергии, J/s = W
    double energyDerivativeW = 0.0;

    // Энергия, соответствующая текущим T, P, Z
    double inventoryEnergyJ = 0.0;

    // Невязка алгебраического уравнения энергии:
    // E_state - inventoryEnergyJ
    double energyResidualJ = 0.0;
};

class WellStirredReactorModel
{
public:
    WellStirredReactorModel(
        WellStirredReactorParameters parameters,
        const ThermoPackage& thermo,
        const ReactionModel& reactions
        );

    const WellStirredReactorParameters& parameters() const;

    WellStirredReactorEvaluation evaluate(
        const ReactorState& state
        ) const;

    EnergyBalanceResult computeEnergyBalance(
        const ReactorState& state
        ) const;

    std::string energyBalanceToString(
        const EnergyBalanceResult& result
        ) const;

    double calculateInventoryEnergyJ(
        const ReactorState& state
        ) const;

    ReactorState initializeEnergyFromEnthalpy(
        const ReactorState& state
        ) const;

    PressureInitializationResult initializePressureForVolume(
        const ReactorState& initialState,
        double pressureMinPa = 1.0e4,
        double pressureMaxPa = 1.0e7,
        double volumeToleranceM3 = 1.0e-8,
        int maxIterations = 100
        ) const;

    PressureInitializationResult initializePressureForVolumeNewton(
        const ReactorState& initialState,
        double pressureMinPa = 1.0e4,
        double pressureMaxPa = 1.0e7,
        double volumeToleranceM3 = 1.0e-8,
        int maxNewtonIterations = 20,
        bool fallbackToBisection = true
        ) const;

    std::string evaluationToString(
        const WellStirredReactorEvaluation& evaluation
        ) const;

    std::string pressureInitializationToString(
        const PressureInitializationResult& result
        ) const;

    MassBalanceResult computeMassBalance(
        const ReactorState& state
        ) const;

    std::string massBalanceToString(
        const MassBalanceResult& result
        ) const;

private:
    WellStirredReactorParameters parameters_;
    const ThermoPackage& thermo_;
    const ReactionModel& reactions_;
};
