#include "mainwindow.h"

#include "core/reactions/ReactionModel.h"
#include "core/ReactorState.h"
#include "core/thermo/ThermoPackage.h"
#include "core/models/WellStirredReactorModel.h"
#include "core/dae/ReducedReactorDae.h"
#include "core/numerics/NewtonSolver.h"

#include <QApplication>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>

void runFlashCase(
    const std::string& title,
    const ThermoPackage& thermo,
    double pressurePa,
    double temperatureK,
    const Composition& z
    )
{
    const MaterialList& materials = thermo.materials();

    FlashResult result = thermo.flash(
        pressurePa,
        temperatureK,
        z
        );

    std::cout << "\n========================================\n";
    std::cout << title << "\n";
    std::cout << "========================================\n";

    std::cout << "P = " << pressurePa << " Pa\n";
    std::cout << "T = " << temperatureK << " K\n";

    std::cout << "\nOverall composition z:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        std::cout << "  z(" << materials[i].name << ") = "
                  << z[i] << "\n";
    }

    std::cout << "\nFlash result:\n";
    std::cout << "  status     = "
              << flashStatusToString(result.status) << "\n";
    std::cout << "  iterations = " << result.iterations << "\n";
    std::cout << "  beta       = " << result.beta << "\n";
    std::cout << "  Z_liq      = " << result.zLiquid << "\n";
    std::cout << "  Z_vap      = " << result.zVapor << "\n";

    std::cout << "\nDensities:\n";
    std::cout << "  rho_mol_liq = "
              << result.molarDensityLiquidKmolPerM3
              << " kmol/m3\n";
    std::cout << "  rho_mol_vap = "
              << result.molarDensityVaporKmolPerM3
              << " kmol/m3\n";
    std::cout << "  rho_mol_mix = "
              << result.molarDensityMixtureKmolPerM3
              << " kmol/m3\n";

    std::cout << "  rho_mass_liq = "
              << result.massDensityLiquidKgPerM3
              << " kg/m3\n";
    std::cout << "  rho_mass_vap = "
              << result.massDensityVaporKgPerM3
              << " kg/m3\n";
    std::cout << "  rho_mass_mix = "
              << result.massDensityMixtureKgPerM3
              << " kg/m3\n";

    std::cout << "\nLiquid composition x:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        std::cout << "  x(" << materials[i].name << ") = "
                  << result.xLiquid[i] << "\n";
    }

    std::cout << "\nVapor composition y:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        std::cout << "  y(" << materials[i].name << ") = "
                  << result.yVapor[i] << "\n";
    }

    std::cout << "\nK-values:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        std::cout << "  K(" << materials[i].name << ") = "
                  << result.kValues[i] << "\n";
    }

    std::cout << "\nMaterial balance check:\n";
    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const double zCalculated =
            (1.0 - result.beta) * result.xLiquid[i]
            + result.beta * result.yVapor[i];

        std::cout << "  z_calc(" << materials[i].name << ") = "
                  << zCalculated
                  << "    error = "
                  << zCalculated - z[i]
                  << "\n";
    }
}

ReactorState explicitEulerMassStep(
    const ReactorState& state,
    const WellStirredReactorModel& reactorModel,
    double timeStepS
    )
{
    if (timeStepS <= 0.0)
    {
        throw std::invalid_argument("Time step must be positive");
    }

    MassBalanceResult massBalance =
        reactorModel.computeMassBalance(state);

    ReactorState newState = state;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        const Component component =
            static_cast<Component>(i);

        const double deltaMass =
            massBalance.massDerivativesKgPerS[i]
            * timeStepS;

        newState.addMassKg(
            component,
            deltaMass
            );
    }

    return newState;
}

void runExplicitEulerMassSimulation(
    const ReactorState& initialState,
    const WellStirredReactorModel& reactorModel,
    double timeStepS,
    int stepCount,
    int printEvery
    )
{
    if (timeStepS <= 0.0)
    {
        throw std::invalid_argument("Time step must be positive");
    }

    if (stepCount <= 0)
    {
        throw std::invalid_argument("Step count must be positive");
    }

    if (printEvery <= 0)
    {
        throw std::invalid_argument("Print interval must be positive");
    }

    ReactorState state = initialState;

    std::cout << "\n========================================\n";
    std::cout << "Explicit Euler mass simulation diagnostic\n";
    std::cout << "========================================\n";

    std::cout << "dt = " << timeStepS << " s\n";
    std::cout << "steps = " << stepCount << "\n";

    std::cout << std::fixed << std::setprecision(6);

    std::cout
        << "\n"
        << "time_s"
        << "\tP_bar"
        << "\tbeta"
        << "\tM_C2H6_kg"
        << "\tM_C5H12_kg"
        << "\tM_H2O_kg"
        << "\tdM_C2H6_dt"
        << "\tdM_C5H12_dt"
        << "\tdM_H2O_dt"
        << "\n";

    for (int step = 0; step <= stepCount; ++step)
    {
        const double timeS =
            step * timeStepS;

        WellStirredReactorEvaluation evaluation =
            reactorModel.evaluate(state);

        MassBalanceResult massBalance =
            reactorModel.computeMassBalance(state);

        if (step % printEvery == 0 || step == stepCount)
        {
            std::cout
                << timeS
                << "\t" << state.pressureBar()
                << "\t" << evaluation.flash.beta
                << "\t" << state.massKg(Component::C2H6)
                << "\t" << state.massKg(Component::C5H12)
                << "\t" << state.massKg(Component::H2O)
                << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::C2H6)]
                << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::C5H12)]
                << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::H2O)]
                << "\n";
        }

        if (step == stepCount)
        {
            break;
        }

        ReactorState predictedState =
            explicitEulerMassStep(
                state,
                reactorModel,
                timeStepS
                );

        PressureInitializationResult pressureInitialization =
            reactorModel.initializePressureForVolumeNewton(
                predictedState,
                1.0e4,
                1.0e7,
                1.0e-6,
                15,
                true
                );

        if (!pressureInitialization.converged)
        {
            std::cout << "\nPressure initialization failed at step "
                      << step + 1 << "\n";
            std::cout << reactorModel.pressureInitializationToString(
                pressureInitialization
                ) << "\n";
            break;
        }

        state =
            pressureInitialization.initializedState;
    }
}

void runImplicitEulerSimulation(
    const ReactorState& initialState,
    const ReducedReactorDae& dae,
    const WellStirredReactorModel& reactorModel,
    double timeStepS,
    int stepCount,
    int printEvery
    )
{
    if (timeStepS <= 0.0)
    {
        throw std::invalid_argument("Time step must be positive");
    }

    if (stepCount <= 0)
    {
        throw std::invalid_argument("Step count must be positive");
    }

    if (printEvery <= 0)
    {
        throw std::invalid_argument("Print interval must be positive");
    }

    ReactorState state = initialState;

    int lastNewtonIterations = 0;

    std::cout << "\n========================================\n";
    std::cout << "Implicit Euler DAE simulation diagnostic\n";
    std::cout << "========================================\n";

    std::cout << "dt = " << timeStepS << " s\n";
    std::cout << "steps = " << stepCount << "\n";

    std::cout << std::fixed << std::setprecision(6);

    std::cout
        << "\n"
        << "time_s"
        << "\tP_bar"
        << "\tbeta"
        << "\tFout"
        << "\tM_C2H6_kg"
        << "\tM_C5H12_kg"
        << "\tM_H2O_kg"
        << "\tdM_C2H6_dt"
        << "\tdM_C5H12_dt"
        << "\tdM_H2O_dt"
        << "\tNewton_it"
        << "\n";

    for (int step = 0; step <= stepCount; ++step)
    {
        const double timeS =
            step * timeStepS;

        WellStirredReactorEvaluation evaluation =
            reactorModel.evaluate(state);

        MassBalanceResult massBalance =
            reactorModel.computeMassBalance(state);

        if (step % printEvery == 0 || step == stepCount)
        {
            std::cout
                << timeS
                << "\t" << state.pressureBar()
                << "\t" << evaluation.flash.beta
                << "\t" << massBalance.outletFlowKmolPerS
                << "\t" << state.massKg(Component::C2H6)
                << "\t" << state.massKg(Component::C5H12)
                << "\t" << state.massKg(Component::H2O)
                << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::C2H6)]
                << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::C5H12)]
                << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::H2O)]
                << "\t" << lastNewtonIterations
                << "\n";
        }

        if (step == stepCount)
        {
            break;
        }

        ReducedReactorDaeStepResult stepResult =
            dae.implicitEulerStep(
                state,
                timeStepS
                );

        if (!stepResult.converged)
        {
            std::cout << "\nImplicit Euler failed at step "
                      << step + 1 << "\n";

            std::cout << dae.stepResultToString(
                stepResult
                ) << "\n";

            break;
        }

        state =
            stepResult.state;

        lastNewtonIterations =
            stepResult.iterations;
    }
}

void runRadauIIA2Simulation(
    const ReactorState& initialState,
    const ReducedReactorDae& dae,
    const WellStirredReactorModel& reactorModel,
    double timeStepS,
    int stepCount,
    int printEvery
    )
{
    if (timeStepS <= 0.0)
    {
        throw std::invalid_argument("Time step must be positive");
    }

    if (stepCount <= 0)
    {
        throw std::invalid_argument("Step count must be positive");
    }

    if (printEvery <= 0)
    {
        throw std::invalid_argument("Print interval must be positive");
    }

    ReactorState state = initialState;

    int lastNewtonIterations = 0;

    std::cout << "\n========================================\n";
    std::cout << "Radau IIA 2-stage DAE simulation diagnostic\n";
    std::cout << "========================================\n";

    std::cout << "dt = " << timeStepS << " s\n";
    std::cout << "steps = " << stepCount << "\n";

    std::cout << std::fixed << std::setprecision(6);

    std::cout
        << "\n"
        << "time_s"
        << "\tP_bar"
        << "\tbeta"
        << "\tFout"
        << "\tM_C2H6_kg"
        << "\tM_C5H12_kg"
        << "\tM_H2O_kg"
        << "\tdM_C2H6_dt"
        << "\tdM_C5H12_dt"
        << "\tdM_H2O_dt"
        << "\tNewton_it"
        << "\n";

    for (int step = 0; step <= stepCount; ++step)
    {
        const double timeS =
            step * timeStepS;

        WellStirredReactorEvaluation evaluation =
            reactorModel.evaluate(state);

        MassBalanceResult massBalance =
            reactorModel.computeMassBalance(state);

        if (step % printEvery == 0 || step == stepCount)
        {
            std::cout
                << timeS
                << "\t" << state.pressureBar()
                << "\t" << evaluation.flash.beta
                << "\t" << massBalance.outletFlowKmolPerS
                << "\t" << state.massKg(Component::C2H6)
                << "\t" << state.massKg(Component::C5H12)
                << "\t" << state.massKg(Component::H2O)
                << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::C2H6)]
                << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::C5H12)]
                << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::H2O)]
                << "\t" << lastNewtonIterations
                << "\n";
        }

        if (step == stepCount)
        {
            break;
        }

        ReducedReactorDaeStepResult stepResult =
            dae.radauIIA2Step(
                state,
                timeStepS
                );

        if (!stepResult.converged)
        {
            std::cout << "\nRadau IIA 2-stage failed at step "
                      << step + 1 << "\n";

            std::cout << dae.stepResultToString(
                stepResult
                ) << "\n";

            break;
        }

        state =
            stepResult.state;

        lastNewtonIterations =
            stepResult.iterations;
    }
}

int main(int argc, char *argv[])
{
    ReactorState state;

    std::cout << state.toString() << std::endl;

    ThermoPackage thermo;

    const MaterialList& materials = thermo.materials();

    std::cout << "\nMaterials:\n";

    for (const Material& material : materials)
    {
        std::cout << "  " << material.name << "\n";
        std::cout << "    Tc    = "
                  << material.criticalTemperatureK << " K\n";
        std::cout << "    Pc    = "
                  << material.criticalPressurePa << " Pa\n";
        std::cout << "    omega = "
                  << material.acentricFactor << "\n";
        std::cout << "    Mw    = "
                  << material.molarMassKgPerKmol << " kg/kmol\n";
    }

    Composition z = state.moleFractions();

    const double pressurePa =
        state.pressureBar() * 1.0e5;

    const double temperatureK =
        state.temperatureC() + 273.15;

    runFlashCase(
        "Case 1: C2H6/C5H12, P = 4 bar, T = 25 C",
        thermo,
        0.4e6,
        25.0 + 273.15,
        Composition{0.5, 0.5, 0.0}
        );

    runFlashCase(
        "Case 2: C2H6/C5H12, P = 20 bar, T = 100 C",
        thermo,
        2.0e6,
        100.0 + 273.15,
        Composition{0.4, 0.6, 0.0}
        );

    runFlashCase(
        "Case 3: reactor initial state",
        thermo,
        pressurePa,
        temperatureK,
        z
        );

    ReactionModel reactionModel(
        thermo.materials()
        );

    // Чтобы строго повторить Aspen-затычку:
    reactionModel.setForwardRateConstant(100.0);
    reactionModel.setReverseRateConstant(100.0);

    // Для мягкой отладки динамики позже можно временно поставить:
    // reactionModel.setForwardRateConstant(0.01);
    // reactionModel.setReverseRateConstant(0.01);

    WellStirredReactorParameters reactorParameters;
    reactorParameters.volumeM3 = 1.0;

    reactorParameters.allowReverseOutletFlow = false;

    WellStirredReactorModel reactorModel(
        reactorParameters,
        thermo,
        reactionModel
        );

    ReducedReactorDae reducedDae(
        reactorModel
        );

    WellStirredReactorEvaluation reactorEvaluation =
        reactorModel.evaluate(state);

    std::cout << "\n========================================\n";
    std::cout << "Reactor model diagnostic\n";
    std::cout << "========================================\n";

    std::cout << reactorModel.evaluationToString(
        reactorEvaluation
        ) << std::endl;

    PressureInitializationResult pressureInitialization =
        reactorModel.initializePressureForVolumeNewton(
            state,
            1.0e4,   // 0.1 bar
            1.0e7,   // 100 bar
            1.0e-8,
            20,
            true     // если Newton не сошелся, использовать бисекцию
            );

    std::cout << "\n========================================\n";
    std::cout << "Pressure initialization\n";
    std::cout << "========================================\n";

    std::cout << reactorModel.pressureInitializationToString(
        pressureInitialization
        ) << std::endl;

    MassBalanceResult massBalance =
        reactorModel.computeMassBalance(
            pressureInitialization.initializedState
            );

    std::cout << "\n========================================\n";
    std::cout << "Mass balance diagnostic\n";
    std::cout << "========================================\n";

    std::cout << reactorModel.massBalanceToString(
        massBalance
        ) << std::endl;

    runImplicitEulerSimulation(
        pressureInitialization.initializedState,
        reducedDae,
        reactorModel,
        1.0,   // dt, s
        30,    // 30 шагов = 30 секунд
        1      // печатать каждый шаг
        );

    const double timeStepS = 1.0;

    ReactorState stateAfterStep =
        explicitEulerMassStep(
            pressureInitialization.initializedState,
            reactorModel,
            timeStepS
            );

    PressureInitializationResult pressureAfterStep =
        reactorModel.initializePressureForVolume(
            stateAfterStep,
            1.0e4,
            1.0e7,
            1.0e-8,
            100
            );

    ReducedReactorDaeVariables daeVariables =
        reducedDae.variablesFromState(
            pressureInitialization.initializedState
            );

    ReducedReactorDaeDerivatives daeDerivatives;

    for (std::size_t i = 0; i < ComponentCount; ++i)
    {
        daeDerivatives.massDerivativesKgPerS[i] =
            massBalance.massDerivativesKgPerS[i];
    }

    ReducedReactorDaeResidual daeResidual =
        reducedDae.computeResidual(
            daeVariables,
            daeDerivatives
            );

    std::cout << "\n========================================\n";
    std::cout << "Reduced DAE residual diagnostic\n";
    std::cout << "========================================\n";

    std::cout << reducedDae.residualToString(
        daeResidual
        ) << std::endl;

    numerics::Vector daeVariableVector =
        reducedDae.variablesToVector(
            daeVariables
            );

    numerics::Vector daeDerivativeVector =
        reducedDae.derivativesToVector(
            daeDerivatives
            );

    numerics::Vector daeResidualVector =
        reducedDae.computeResidualVector(
            daeVariableVector,
            daeDerivativeVector,
            daeVariables.temperatureC
            );

    std::cout << "\n========================================\n";
    std::cout << "Reduced DAE vector residual diagnostic\n";
    std::cout << "========================================\n";

    std::cout << std::fixed << std::setprecision(10);

    for (std::size_t i = 0; i < daeResidualVector.size(); ++i)
    {
        std::cout << "R[" << i << "] = "
                  << daeResidualVector[i]
                  << "\n";
    }

    ReducedReactorDaeStepResult implicitStep =
        reducedDae.implicitEulerStep(
            pressureInitialization.initializedState,
            1.0
            );

    std::cout << "\n========================================\n";
    std::cout << "Reduced DAE implicit Euler diagnostic\n";
    std::cout << "========================================\n";

    std::cout << reducedDae.stepResultToString(
        implicitStep
        ) << std::endl;

    MassBalanceResult implicitStepMassBalance =
        reactorModel.computeMassBalance(
            implicitStep.state
            );

    std::cout << "\nMass balance after implicit Euler step:\n";
    std::cout << reactorModel.massBalanceToString(
        implicitStepMassBalance
        ) << std::endl;

    ReducedReactorDaeStepResult radauStep =
        reducedDae.radauIIA2Step(
            pressureInitialization.initializedState,
            1.0
            );

    std::cout << "\n========================================\n";
    std::cout << "Reduced DAE Radau IIA 2-stage diagnostic\n";
    std::cout << "========================================\n";

    std::cout << reducedDae.stepResultToString(
        radauStep
        ) << std::endl;

    MassBalanceResult radauStepMassBalance =
        reactorModel.computeMassBalance(
            radauStep.state
            );

    std::cout << "\nMass balance after Radau IIA 2-stage step:\n";
    std::cout << reactorModel.massBalanceToString(
        radauStepMassBalance
        ) << std::endl;

    runRadauIIA2Simulation(
        pressureInitialization.initializedState,
        reducedDae,
        reactorModel,
        1.0,   // dt, s
        30,    // 30 шагов = 30 секунд
        1      // печатать каждый шаг
        );


    {
        numerics::NewtonOptions options;
        options.residualTolerance = 1.0e-12;
        options.stepTolerance = 1.0e-12;
        options.maxIterations = 30;

        numerics::Vector initialGuess{1.0};

        auto function = [](const numerics::Vector& x)
        {
            numerics::Vector f(1);

            // Решаем x^2 - 2 = 0.
            // Правильный положительный корень sqrt(2) ≈ 1.41421356.
            f[0] = x[0] * x[0] - 2.0;

            return f;
        };

        numerics::NewtonResult result =
            numerics::NewtonSolver::solve(
                function,
                initialGuess,
                options
                );

        std::cout << "\n========================================\n";
        std::cout << "Newton solver simple diagnostic\n";
        std::cout << "========================================\n";

        std::cout << "converged     = "
                  << (result.converged ? "true" : "false") << "\n";
        std::cout << "iterations    = "
                  << result.iterations << "\n";
        std::cout << "solution x    = "
                  << result.solution[0] << "\n";
        std::cout << "residual      = "
                  << result.residual[0] << "\n";
        std::cout << "residual norm = "
                  << result.residualNorm << "\n";
        std::cout << "message       = "
                  << result.message << "\n";
    }


    std::cout << "\n========================================\n";
    std::cout << "One explicit Euler mass step diagnostic\n";
    std::cout << "========================================\n";

    std::cout << "time step = "
              << timeStepS
              << " s\n";

    std::cout << reactorModel.pressureInitializationToString(
        pressureAfterStep
        ) << std::endl;

    MassBalanceResult massBalanceAfterStep =
        reactorModel.computeMassBalance(
            pressureAfterStep.initializedState
            );

    std::cout << "\nMass balance after one step:\n";
    std::cout << reactorModel.massBalanceToString(
        massBalanceAfterStep
        ) << std::endl;

    runExplicitEulerMassSimulation(
        pressureInitialization.initializedState,
        reactorModel,
        0.005,  // dt, s
        200,    // 200 шагов = 1 секунда
        20      // печатать каждую 20-ю строку
        );

    QApplication app(argc, argv);

    MainWindow window;
    window.show();

    return app.exec();
}
