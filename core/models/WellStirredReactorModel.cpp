#include "WellStirredReactorModel.h"
#include "../numerics/NewtonSolver.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>

WellStirredReactorModel::WellStirredReactorModel(
    WellStirredReactorParameters parameters,
    const ThermoPackage& thermo,
    const ReactionModel& reactions
    )
    : parameters_(parameters),
    thermo_(thermo),
    reactions_(reactions)
{
    if (parameters_.volumeM3 <= 0.0)
    {
        throw std::invalid_argument("Reactor volume must be positive");
    }
}

const WellStirredReactorParameters&
WellStirredReactorModel::parameters() const
{
    return parameters_;
}

WellStirredReactorEvaluation WellStirredReactorModel::evaluate(
    const ReactorState& state
    ) const
{
    WellStirredReactorEvaluation evaluation;

    evaluation.pressurePa =
        state.pressureBar() * 1.0e5;

    evaluation.temperatureK =
        state.temperatureC() + 273.15;

    evaluation.totalMolesKmol =
        state.totalMolesKmol();

    evaluation.zOverall =
        state.moleFractions();

    evaluation.flash =
        thermo_.flash(
            evaluation.pressurePa,
            evaluation.temperatureK,
            evaluation.zOverall
            );

    if (evaluation.flash.molarDensityMixtureKmolPerM3 <= 0.0)
    {
        throw std::runtime_error("Mixture molar density must be positive");
    }

    evaluation.calculatedVolumeM3 =
        evaluation.totalMolesKmol
        / evaluation.flash.molarDensityMixtureKmolPerM3;

    evaluation.volumeResidualM3 =
        evaluation.calculatedVolumeM3
        - parameters_.volumeM3;

    return evaluation;
}

double WellStirredReactorModel::calculateInventoryEnergyJ(
    const ReactorState& state
    ) const
{
    const WellStirredReactorEvaluation evaluation =
        evaluate(state);

    return
        evaluation.totalMolesKmol
        * evaluation.flash.molarEnthalpyMixtureJPerKmol;
}



ReactorState WellStirredReactorModel::initializeEnergyFromEnthalpy(
    const ReactorState& state
    ) const
{
    ReactorState initializedState =
        state;

    const double energy =
        calculateInventoryEnergyJ(state);

    initializedState.setEnergyJ(
        energy
        );

    return initializedState;
}

PressureInitializationResult
WellStirredReactorModel::initializePressureForVolume(
    const ReactorState& initialState,
    double pressureMinPa,
    double pressureMaxPa,
    double volumeToleranceM3,
    int maxIterations
    ) const
{
    if (pressureMinPa <= 0.0)
    {
        throw std::invalid_argument("Minimum pressure must be positive");
    }

    if (pressureMaxPa <= pressureMinPa)
    {
        throw std::invalid_argument("Maximum pressure must be greater than minimum pressure");
    }

    if (volumeToleranceM3 <= 0.0)
    {
        throw std::invalid_argument("Volume tolerance must be positive");
    }

    if (maxIterations <= 0)
    {
        throw std::invalid_argument("Maximum iteration count must be positive");
    }

    auto evaluateAtPressure =
        [this, &initialState](double pressurePa)
    {
        ReactorState state = initialState;
        state.setPressureBar(pressurePa / 1.0e5);
        return evaluate(state);
    };

    PressureInitializationResult result;
    result.initializedState = initialState;

    WellStirredReactorEvaluation evaluationLow =
        evaluateAtPressure(pressureMinPa);

    WellStirredReactorEvaluation evaluationHigh =
        evaluateAtPressure(pressureMaxPa);

    double residualLow = evaluationLow.volumeResidualM3;
    double residualHigh = evaluationHigh.volumeResidualM3;

    if (std::abs(residualLow) < volumeToleranceM3)
    {
        result.converged = true;
        result.iterations = 0;
        result.pressurePa = pressureMinPa;
        result.pressureBar = pressureMinPa / 1.0e5;
        result.residualM3 = residualLow;
        result.initializedState.setPressureBar(result.pressureBar);
        result.evaluation = evaluationLow;
        return result;
    }

    if (std::abs(residualHigh) < volumeToleranceM3)
    {
        result.converged = true;
        result.iterations = 0;
        result.pressurePa = pressureMaxPa;
        result.pressureBar = pressureMaxPa / 1.0e5;
        result.residualM3 = residualHigh;
        result.initializedState.setPressureBar(result.pressureBar);
        result.evaluation = evaluationHigh;
        return result;
    }

    const bool sameSign =
        (residualLow > 0.0 && residualHigh > 0.0)
        || (residualLow < 0.0 && residualHigh < 0.0);

    if (sameSign)
    {
        // Не удалось зажать корень на заданном интервале.
        // Возвращаем диагностику для исходного давления.

        result.converged = false;
        result.iterations = 0;
        result.pressurePa = initialState.pressureBar() * 1.0e5;
        result.pressureBar = initialState.pressureBar();
        result.evaluation = evaluate(initialState);
        result.residualM3 = result.evaluation.volumeResidualM3;
        return result;
    }

    double leftPressure = pressureMinPa;
    double rightPressure = pressureMaxPa;

    double leftResidual = residualLow;

    for (int iteration = 1; iteration <= maxIterations; ++iteration)
    {
        const double middlePressure =
            0.5 * (leftPressure + rightPressure);

        WellStirredReactorEvaluation middleEvaluation =
            evaluateAtPressure(middlePressure);

        const double middleResidual =
            middleEvaluation.volumeResidualM3;



        if (std::abs(middleResidual) < volumeToleranceM3)
        {
            result.converged = true;
            result.iterations = iteration;
            result.pressurePa = middlePressure;
            result.pressureBar = middlePressure / 1.0e5;
            result.residualM3 = middleResidual;
            result.initializedState = initialState;
            result.initializedState.setPressureBar(result.pressureBar);
            result.evaluation = middleEvaluation;
            return result;
        }

        const bool rootBetweenLeftAndMiddle =
            (leftResidual > 0.0 && middleResidual < 0.0)
            || (leftResidual < 0.0 && middleResidual > 0.0);

        if (rootBetweenLeftAndMiddle)
        {
            rightPressure = middlePressure;
        }
        else
        {
            leftPressure = middlePressure;
            leftResidual = middleResidual;
        }
    }

    const double finalPressure =
        0.5 * (leftPressure + rightPressure);

    WellStirredReactorEvaluation finalEvaluation =
        evaluateAtPressure(finalPressure);

    result.converged = false;
    result.iterations = maxIterations;
    result.pressurePa = finalPressure;
    result.pressureBar = finalPressure / 1.0e5;
    result.residualM3 = finalEvaluation.volumeResidualM3;
    result.initializedState = initialState;
    result.initializedState.setPressureBar(result.pressureBar);
    result.evaluation = finalEvaluation;

    return result;
}

PressureInitializationResult
WellStirredReactorModel::initializePressureForVolumeNewton(
    const ReactorState& initialState,
    double pressureMinPa,
    double pressureMaxPa,
    double volumeToleranceM3,
    int maxNewtonIterations,
    bool fallbackToBisection
    ) const
{
    if (pressureMinPa <= 0.0)
    {
        throw std::invalid_argument("Minimum pressure must be positive");
    }

    if (pressureMaxPa <= pressureMinPa)
    {
        throw std::invalid_argument("Maximum pressure must be greater than minimum pressure");
    }

    if (volumeToleranceM3 <= 0.0)
    {
        throw std::invalid_argument("Volume tolerance must be positive");
    }

    if (maxNewtonIterations <= 0)
    {
        throw std::invalid_argument("Maximum Newton iteration count must be positive");
    }

    const double pressureMinBar =
        pressureMinPa / 1.0e5;

    const double pressureMaxBar =
        pressureMaxPa / 1.0e5;

    double initialPressureBar =
        initialState.pressureBar();

    initialPressureBar =
        std::max(
            pressureMinBar,
            std::min(pressureMaxBar, initialPressureBar)
            );

    auto evaluateAtPressureBar =
        [this, &initialState](double pressureBar)
    {
        ReactorState state = initialState;
        state.setPressureBar(pressureBar);
        return evaluate(state);
    };

    auto residualFunction =
        [evaluateAtPressureBar](const numerics::Vector& x)
    {
        if (x.size() != 1)
        {
            throw std::runtime_error("Pressure Newton function expects one variable");
        }

        const double pressureBar =
            x[0];

        WellStirredReactorEvaluation evaluation =
            evaluateAtPressureBar(pressureBar);

        return numerics::Vector{
            evaluation.volumeResidualM3
        };
    };

    try
    {
        numerics::NewtonOptions options;

        options.maxIterations =
            maxNewtonIterations;

        options.residualTolerance =
            volumeToleranceM3;

        options.stepTolerance =
            1.0e-10;

        options.finiteDifferenceStep =
            1.0e-6;

        options.minDampingFactor =
            1.0e-6;

        options.lowerBounds =
            numerics::Vector{pressureMinBar};

        options.upperBounds =
            numerics::Vector{pressureMaxBar};

        numerics::NewtonResult newtonResult =
            numerics::NewtonSolver::solve(
                residualFunction,
                numerics::Vector{initialPressureBar},
                options
                );

        if (newtonResult.converged)
        {
            const double pressureBar =
                newtonResult.solution[0];

            WellStirredReactorEvaluation finalEvaluation =
                evaluateAtPressureBar(pressureBar);

            PressureInitializationResult result;

            result.converged = true;
            result.iterations = newtonResult.iterations;
            result.pressureBar = pressureBar;
            result.pressurePa = pressureBar * 1.0e5;
            result.residualM3 = finalEvaluation.volumeResidualM3;
            result.method = "newton";
            result.message = newtonResult.message;

            result.initializedState = initialState;
            result.initializedState.setPressureBar(pressureBar);

            result.evaluation = finalEvaluation;

            return result;
        }

        if (!fallbackToBisection)
        {
            PressureInitializationResult result;

            result.converged = false;
            result.iterations = newtonResult.iterations;
            result.pressureBar = newtonResult.solution.empty()
                                     ? initialPressureBar
                                     : newtonResult.solution[0];
            result.pressurePa = result.pressureBar * 1.0e5;
            result.method = "newton";
            result.message = newtonResult.message;

            result.initializedState = initialState;
            result.initializedState.setPressureBar(result.pressureBar);

            result.evaluation =
                evaluateAtPressureBar(result.pressureBar);

            result.residualM3 =
                result.evaluation.volumeResidualM3;

            return result;
        }

        PressureInitializationResult fallbackResult =
            initializePressureForVolume(
                initialState,
                pressureMinPa,
                pressureMaxPa,
                volumeToleranceM3,
                100
                );

        fallbackResult.method =
            "newton_failed_then_bisection";

        fallbackResult.message =
            "Newton failed: " + newtonResult.message;

        return fallbackResult;
    }
    catch (const std::exception& exception)
    {
        if (!fallbackToBisection)
        {
            throw;
        }

        PressureInitializationResult fallbackResult =
            initializePressureForVolume(
                initialState,
                pressureMinPa,
                pressureMaxPa,
                volumeToleranceM3,
                100
                );

        fallbackResult.method =
            "newton_exception_then_bisection";

        fallbackResult.message =
            exception.what();

        return fallbackResult;
    }
}

MassBalanceResult WellStirredReactorModel::computeMassBalance(
    const ReactorState& state
    ) const
{
    MassBalanceResult result;

    result.evaluation = evaluate(state);

    const MaterialList& materials =
        thermo_.materials();

    result.inletFlowKmolPerS =
        parameters_.inletFlowKmolPerS;

    result.inletComposition =
        parameters_.inletComposition;

    // Как в Aspen:
    // Outlet.Z = y
    result.outletComposition =
        result.evaluation.flash.yVapor;

    // Пока повторяем Aspen-затычку:
    // Outlet.F = 0.1 * (P - 2.5)
    result.rawOutletFlowKmolPerS =
        parameters_.outletValveCoefficientKmolPerSBar
        * (
            state.pressureBar()
            - parameters_.outletReferencePressureBar
            );

    if (parameters_.allowReverseOutletFlow)
    {
        result.outletFlowKmolPerS =
            result.rawOutletFlowKmolPerS;
    }
    else
    {
        result.outletFlowKmolPerS =
            std::max(
                0.0,
                result.rawOutletFlowKmolPerS
                );
    }

    const double beta =
        result.evaluation.flash.beta;

    const double rhoLiquid =
        result.evaluation.flash.molarDensityLiquidKmolPerM3;

    const double rhoMix =
        result.evaluation.flash.molarDensityMixtureKmolPerM3;

    if (rhoLiquid > 0.0)
    {
        result.liquidVolumeFraction =
            (1.0 - beta)
            / rhoLiquid
            * rhoMix;
    }
    else
    {
        result.liquidVolumeFraction = 0.0;
    }

    result.reactionRates =
        reactions_.computeRates(
            result.evaluation.temperatureK,
            result.evaluation.pressurePa,
            result.evaluation.flash.xLiquid,
            result.evaluation.flash.molarDensityLiquidKmolPerM3
            );

    result.totalInletMassFlowKgPerS = 0.0;
    result.totalOutletMassFlowKgPerS = 0.0;
    result.totalReactionMassRateKgPerS = 0.0;
    result.totalMassDerivativeKgPerS = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const double molecularWeight =
            materials[i].molarMassKgPerKmol;

        result.inletMassFlowKgPerS[i] =
            result.inletFlowKmolPerS
            * result.inletComposition[i]
            * molecularWeight;

        result.outletMassFlowKgPerS[i] =
            result.outletFlowKmolPerS
            * result.outletComposition[i]
            * molecularWeight;

        // Аналог Aspen:
        //
        // Reactions.Rate * Mw * V * Fliq
        //
        // componentRatesKmolPerM3S: kmol/(m3_liquid*s)
        // molecularWeight: kg/kmol
        // volumeM3 * liquidVolumeFraction: m3_liquid
        //
        // результат: kg/s
        result.reactionMassRateKgPerS[i] =
            result.reactionRates.componentRatesKmolPerM3S[i]
            * molecularWeight
            * parameters_.volumeM3
            * result.liquidVolumeFraction;

        result.massDerivativesKgPerS[i] =
            result.inletMassFlowKgPerS[i]
            - result.outletMassFlowKgPerS[i]
            + result.reactionMassRateKgPerS[i];

        result.totalInletMassFlowKgPerS +=
            result.inletMassFlowKgPerS[i];

        result.totalOutletMassFlowKgPerS +=
            result.outletMassFlowKgPerS[i];

        result.totalReactionMassRateKgPerS +=
            result.reactionMassRateKgPerS[i];

        result.totalMassDerivativeKgPerS +=
            result.massDerivativesKgPerS[i];
    }

    return result;
}

EnergyBalanceResult WellStirredReactorModel::computeEnergyBalance(
    const ReactorState& state
    ) const
{
    EnergyBalanceResult result;

    result.evaluation =
        evaluate(state);

    const double reactorTemperatureK =
        result.evaluation.temperatureK;

    const double inletPressurePa =
        parameters_.inletPressureBar * 1.0e5;

    const double inletTemperatureK =
        parameters_.inletTemperatureC + 273.15;

    const FlashResult inletFlash =
        thermo_.flash(
            inletPressurePa,
            inletTemperatureK,
            parameters_.inletComposition
            );

    result.inletMolarEnthalpyJPerKmol =
        inletFlash.molarEnthalpyMixtureJPerKmol;

    result.reactorMolarEnthalpyJPerKmol =
        result.evaluation.flash.molarEnthalpyMixtureJPerKmol;

    result.outletMolarEnthalpyJPerKmol =
        result.evaluation.flash.molarEnthalpyVaporJPerKmol;

    const MassBalanceResult massBalance =
        computeMassBalance(state);

    result.inletEnthalpyFlowW =
        parameters_.inletFlowKmolPerS
        * result.inletMolarEnthalpyJPerKmol;

    result.outletEnthalpyFlowW =
        massBalance.outletFlowKmolPerS
        * result.outletMolarEnthalpyJPerKmol;

    const double reactorTemperatureC =
        state.temperatureC();

    result.heatTransferW =
        parameters_.heatTransferCoefficientWPerM2K
        * parameters_.heatTransferAreaM2
        * (
            parameters_.jacketTemperatureC
            - reactorTemperatureC
            );

    result.energyDerivativeW =
        result.inletEnthalpyFlowW
        - result.outletEnthalpyFlowW
        + result.heatTransferW;

    result.inventoryEnergyJ =
        result.evaluation.totalMolesKmol
        * result.reactorMolarEnthalpyJPerKmol;

    result.energyResidualJ =
        state.energyJ()
        - result.inventoryEnergyJ;

    return result;
}

std::string WellStirredReactorModel::evaluationToString(
    const WellStirredReactorEvaluation& evaluation
    ) const
{
    const MaterialList& materials =
        thermo_.materials();

    std::ostringstream out;

    out << std::fixed << std::setprecision(8);

    out << "Well-stirred reactor evaluation:\n";

    out << "  Conditions:\n";
    out << "    P = " << evaluation.pressurePa << " Pa\n";
    out << "    T = " << evaluation.temperatureK << " K\n";

    out << "  Total amount:\n";
    out << "    N = " << evaluation.totalMolesKmol << " kmol\n";

    out << "  Overall composition Z:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "    Z(" << materials[i].name << ") = "
            << evaluation.zOverall[i] << "\n";
    }

    out << "  Flash:\n";
    out << "    status = "
        << flashStatusToString(evaluation.flash.status) << "\n";
    out << "    beta   = "
        << evaluation.flash.beta << "\n";
    out << "    Z_liq  = "
        << evaluation.flash.zLiquid << "\n";
    out << "    Z_vap  = "
        << evaluation.flash.zVapor << "\n";

    out << "  Phase compositions:\n";

    out << "    Liquid x:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "      x(" << materials[i].name << ") = "
            << evaluation.flash.xLiquid[i] << "\n";
    }

    out << "    Vapor y:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "      y(" << materials[i].name << ") = "
            << evaluation.flash.yVapor[i] << "\n";
    }

    out << "  Densities:\n";
    out << "    rho_mol_liq  = "
        << evaluation.flash.molarDensityLiquidKmolPerM3
        << " kmol/m3\n";
    out << "    rho_mol_vap  = "
        << evaluation.flash.molarDensityVaporKmolPerM3
        << " kmol/m3\n";
    out << "    rho_mol_mix  = "
        << evaluation.flash.molarDensityMixtureKmolPerM3
        << " kmol/m3\n";

    out << "    rho_mass_liq = "
        << evaluation.flash.massDensityLiquidKgPerM3
        << " kg/m3\n";
    out << "    rho_mass_vap = "
        << evaluation.flash.massDensityVaporKgPerM3
        << " kg/m3\n";
    out << "    rho_mass_mix = "
        << evaluation.flash.massDensityMixtureKgPerM3
        << " kg/m3\n";

    out << "  Volume equation:\n";
    out << "    specified V  = "
        << parameters_.volumeM3 << " m3\n";
    out << "    calculated V = "
        << evaluation.calculatedVolumeM3 << " m3\n";
    out << "    residual     = "
        << evaluation.volumeResidualM3 << " m3\n";

    return out.str();
}

std::string WellStirredReactorModel::pressureInitializationToString(
    const PressureInitializationResult& result
    ) const
{
    std::ostringstream out;

    out << std::fixed << std::setprecision(8);

    out << "Pressure initialization result:\n";
    out << "  converged  = "
        << (result.converged ? "true" : "false") << "\n";
    out << "  iterations = "
        << result.iterations << "\n";
    out << "  method     = "
        << result.method << "\n";
    out << "  message    = "
        << result.message << "\n";
    out << "  pressure   = "
        << result.pressurePa << " Pa\n";
    out << "  pressure   = "
        << result.pressureBar << " bar\n";
    out << "  residual   = "
        << result.residualM3 << " m3\n";

    out << "\nInitialized reactor state:\n";
    out << result.initializedState.toString();

    out << "\nEvaluation at initialized pressure:\n";
    out << evaluationToString(result.evaluation);

    return out.str();
}

std::string WellStirredReactorModel::massBalanceToString(
    const MassBalanceResult& result
    ) const
{
    const MaterialList& materials =
        thermo_.materials();

    std::ostringstream out;

    out << std::fixed << std::setprecision(8);

    out << "Mass balance result:\n";

    out << "  Flows:\n";
    out << "    inlet F      = "
        << result.inletFlowKmolPerS
        << " kmol/s\n";
    out << "    raw outlet F = "
        << result.rawOutletFlowKmolPerS
        << " kmol/s\n";
    out << "    outlet F     = "
        << result.outletFlowKmolPerS
        << " kmol/s\n";

    out << "  Liquid volume fraction:\n";
    out << "    Fliq = "
        << result.liquidVolumeFraction
        << "\n";

    out << "  Inlet composition:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "    Zin(" << materials[i].name << ") = "
            << result.inletComposition[i] << "\n";
    }

    out << "  Outlet composition:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "    Zout(" << materials[i].name << ") = "
            << result.outletComposition[i] << "\n";
    }

    out << "  Reaction diagnostics:\n";
    out << "    r_forward = "
        << result.reactionRates.forwardRateKmolPerM3S
        << " kmol/(m3_liq*s)\n";
    out << "    r_reverse = "
        << result.reactionRates.reverseRateKmolPerM3S
        << " kmol/(m3_liq*s)\n";

    out << "  Component source terms from reaction:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "    R(" << materials[i].name << ") = "
            << result.reactionRates.componentRatesKmolPerM3S[i]
            << " kmol/(m3_liq*s)\n";
    }

    out << "  Mass flow contributions:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        out << "    " << materials[i].name << ":\n";
        out << "      inlet    = "
            << result.inletMassFlowKgPerS[i]
            << " kg/s\n";
        out << "      outlet   = "
            << result.outletMassFlowKgPerS[i]
            << " kg/s\n";
        out << "      reaction = "
            << result.reactionMassRateKgPerS[i]
            << " kg/s\n";
        out << "      dM/dt    = "
            << result.massDerivativesKgPerS[i]
            << " kg/s\n";
    }

    out << "  Total mass contributions:\n";
    out << "    inlet total    = "
        << result.totalInletMassFlowKgPerS
        << " kg/s\n";
    out << "    outlet total   = "
        << result.totalOutletMassFlowKgPerS
        << " kg/s\n";
    out << "    reaction total = "
        << result.totalReactionMassRateKgPerS
        << " kg/s\n";
    out << "    dM(total)/dt   = "
        << result.totalMassDerivativeKgPerS
        << " kg/s\n";

    return out.str();
}

std::string WellStirredReactorModel::energyBalanceToString(
    const EnergyBalanceResult& result
    ) const
{
    std::ostringstream out;

    out << std::fixed << std::setprecision(8);

    out << "Energy balance result:\n";

    out << "  Molar enthalpies:\n";
    out << "    h_inlet   = "
        << result.inletMolarEnthalpyJPerKmol
        << " J/kmol\n";
    out << "    h_outlet  = "
        << result.outletMolarEnthalpyJPerKmol
        << " J/kmol\n";
    out << "    h_reactor = "
        << result.reactorMolarEnthalpyJPerKmol
        << " J/kmol\n";

    out << "  Energy flows:\n";
    out << "    Hdot_in   = "
        << result.inletEnthalpyFlowW
        << " W\n";
    out << "    Hdot_out  = "
        << result.outletEnthalpyFlowW
        << " W\n";
    out << "    Qdot      = "
        << result.heatTransferW
        << " W\n";
    out << "    dE/dt     = "
        << result.energyDerivativeW
        << " W\n";

    out << "  Inventory energy:\n";
    out << "    E_state   = "
        << result.evaluation.totalMolesKmol
        << " kmol * h_mix = "
        << result.inventoryEnergyJ
        << " J\n";
    out << "    residual  = "
        << result.energyResidualJ
        << " J\n";

    return out.str();
}
