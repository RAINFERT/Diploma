#include "mainwindow.h"

#include "core/ReactorState.h"
#include "core/dae/ReducedReactorDae.h"
#include "core/models/WellStirredReactorModel.h"
#include "core/numerics/NewtonSolver.h"
#include "core/reactions/ReactionModel.h"
#include "core/thermo/ThermoPackage.h"
#include "core/simulation/ReactorSimulation.h"
#include "core/dae/EnergyReactorDae.h"
#include "core/simulation/EnergyReactorSimulation.h"

#include <QApplication>

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <cmath>

namespace
{
    // ============================================================
    // Output switches
    // ============================================================
    // Keep diagnostic functions in the code, but disable most of them
    // during a normal run. Turn a flag to true when you need to check
    // a specific subsystem again.

    constexpr bool PrintInitialState = true;
    constexpr bool PrintMaterials = false;

    constexpr bool RunFlashDiagnostics = true;
    constexpr bool RunReactorModelDiagnostic = false;
    constexpr bool RunPressureInitializationDiagnostic = true;
    constexpr bool RunMassBalanceDiagnostic = false;

    constexpr bool RunDaeResidualDiagnostic = false;
    constexpr bool RunNewtonSimpleDiagnostic = false;

    constexpr bool RunImplicitEulerSingleStepDiagnostic = false;
    constexpr bool RunImplicitEulerSimulationDiagnostic = false;

    constexpr bool RunRadauSingleStepDiagnostic = false;
    constexpr bool RunRadauSimulationDiagnostic = false;

    constexpr bool RunClosedHeatingValidation = false;

    constexpr bool RunExplicitEulerSingleStepDiagnostic = false;
    constexpr bool RunExplicitEulerSimulationDiagnostic = false;

    constexpr bool RunEnthalpyDiagnostic = false;
    constexpr bool RunEnergyBalanceDiagnostic = false;
    constexpr bool RunEnergyDaeResidualDiagnostic = false;
    constexpr bool RunEnergyRadauSingleStepDiagnostic = false;
    constexpr bool RunEnergyRadauSimulationDiagnostic = false;

    constexpr bool PrintFinalState = true;

    // ============================================================
    // Simulation settings
    // ============================================================

    constexpr double InitialPressureSearchMinPa = 1.0e4;  // 0.1 bar
    constexpr double InitialPressureSearchMaxPa = 1.0e7;  // 100 bar
    constexpr double InitialVolumeToleranceM3 = 1.0e-8;
    constexpr int InitialPressureNewtonMaxIterations = 20;

    constexpr double MainTimeStepS = 1.0;
    constexpr int MainStepCount = 30;
    constexpr int MainPrintEvery = 1;

    // ============================================================
    // Small printing helpers
    // ============================================================

    void printSectionTitle(const std::string& title)
    {
        std::cout << "\n========================================\n";
        std::cout << title << "\n";
        std::cout << "========================================\n";
    }

    void printMaterials(const ThermoPackage& thermo)
    {
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
    }

    void printPressureInitializationSummary(
        const PressureInitializationResult& result
    )
    {
        printSectionTitle("Pressure initialization summary");

        std::cout << std::fixed << std::setprecision(8);
        std::cout << "converged  = "
                  << (result.converged ? "true" : "false") << "\n";
        std::cout << "method     = " << result.method << "\n";
        std::cout << "iterations = " << result.iterations << "\n";
        std::cout << "P          = " << result.pressureBar << " bar\n";
        std::cout << "residual   = " << result.residualM3 << " m3\n";
        std::cout << "message    = " << result.message << "\n";
    }

    // ============================================================
    // Flash diagnostics
    // ============================================================

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

        printSectionTitle(title);

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

        std::cout << "\nEnthalpies:\n";
        std::cout << "  H_liq = "
                  << result.molarEnthalpyLiquidJPerKmol
                  << " J/kmol\n";
        std::cout << "  H_vap = "
                  << result.molarEnthalpyVaporJPerKmol
                  << " J/kmol\n";
        std::cout << "  H_mix = "
                  << result.molarEnthalpyMixtureJPerKmol
                  << " J/kmol\n";

        std::cout << "  H_liq = "
                  << result.molarEnthalpyLiquidJPerKmol / 1000.0
                  << " kJ/kmol\n";
        std::cout << "  H_vap = "
                  << result.molarEnthalpyVaporJPerKmol / 1000.0
                  << " kJ/kmol\n";
        std::cout << "  H_mix = "
                  << result.molarEnthalpyMixtureJPerKmol / 1000.0
                  << " kJ/kmol\n";

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

    void runFlashDiagnostics(
        const ReactorState& state,
        const ThermoPackage& thermo
    )
    {
        const Composition z = state.moleFractions();

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
    }

    // ============================================================
    // Explicit Euler diagnostics
    // ============================================================

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

        const MassBalanceResult massBalance =
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

    ReactorState runExplicitEulerMassSimulation(
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

        printSectionTitle("Explicit Euler mass simulation diagnostic");

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

            const WellStirredReactorEvaluation evaluation =
                reactorModel.evaluate(state);

            const MassBalanceResult massBalance =
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

            const PressureInitializationResult pressureInitialization =
                reactorModel.initializePressureForVolumeNewton(
                    predictedState,
                    InitialPressureSearchMinPa,
                    InitialPressureSearchMaxPa,
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

        return state;
    }

    // ============================================================
    // Implicit Euler diagnostics
    // ============================================================

    ReactorState runImplicitEulerSimulation(
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

        printSectionTitle("Implicit Euler DAE simulation diagnostic");

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

            const WellStirredReactorEvaluation evaluation =
                reactorModel.evaluate(state);

            const MassBalanceResult massBalance =
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

            const ReducedReactorDaeStepResult stepResult =
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

            state = stepResult.state;
            lastNewtonIterations = stepResult.iterations;
        }

        return state;
    }

    // ============================================================
    // Isolated diagnostics
    // ============================================================

    void runNewtonSimpleDiagnostic()
    {
        numerics::NewtonOptions options;
        options.residualTolerance = 1.0e-12;
        options.stepTolerance = 1.0e-12;
        options.maxIterations = 30;

        const numerics::Vector initialGuess{1.0};

        auto function = [](const numerics::Vector& x)
        {
            numerics::Vector f(1);
            f[0] = x[0] * x[0] - 2.0;
            return f;
        };

        const numerics::NewtonResult result =
            numerics::NewtonSolver::solve(
                function,
                initialGuess,
                options
            );

        printSectionTitle("Newton solver simple diagnostic");

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

    void runDaeResidualDiagnostics(
        const ReactorState& initializedState,
        const ReducedReactorDae& reducedDae,
        const MassBalanceResult& massBalance
    )
    {
        ReducedReactorDaeVariables daeVariables =
            reducedDae.variablesFromState(initializedState);

        ReducedReactorDaeDerivatives daeDerivatives;

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            daeDerivatives.massDerivativesKgPerS[i] =
                massBalance.massDerivativesKgPerS[i];
        }

        const ReducedReactorDaeResidual daeResidual =
            reducedDae.computeResidual(
                daeVariables,
                daeDerivatives
            );

        printSectionTitle("Reduced DAE residual diagnostic");
        std::cout << reducedDae.residualToString(daeResidual) << std::endl;

        const numerics::Vector daeVariableVector =
            reducedDae.variablesToVector(daeVariables);

        const numerics::Vector daeDerivativeVector =
            reducedDae.derivativesToVector(daeDerivatives);

        const numerics::Vector daeResidualVector =
            reducedDae.computeResidualVector(
                daeVariableVector,
                daeDerivativeVector,
                daeVariables.temperatureC
            );

        printSectionTitle("Reduced DAE vector residual diagnostic");

        std::cout << std::fixed << std::setprecision(10);

        for (std::size_t i = 0; i < daeResidualVector.size(); ++i)
        {
            std::cout << "R[" << i << "] = "
                      << daeResidualVector[i]
                      << "\n";
        }
    }

    void runExplicitEulerSingleStepDiagnostic(
        const ReactorState& initializedState,
        const WellStirredReactorModel& reactorModel
    )
    {
        constexpr double timeStepS = 1.0;

        const ReactorState stateAfterStep =
            explicitEulerMassStep(
                initializedState,
                reactorModel,
                timeStepS
            );

        const PressureInitializationResult pressureAfterStep =
            reactorModel.initializePressureForVolumeNewton(
                stateAfterStep,
                InitialPressureSearchMinPa,
                InitialPressureSearchMaxPa,
                InitialVolumeToleranceM3,
                InitialPressureNewtonMaxIterations,
                true
            );

        printSectionTitle("One explicit Euler mass step diagnostic");

        std::cout << "time step = "
                  << timeStepS
                  << " s\n";

        std::cout << reactorModel.pressureInitializationToString(
            pressureAfterStep
        ) << std::endl;

        const MassBalanceResult massBalanceAfterStep =
            reactorModel.computeMassBalance(
                pressureAfterStep.initializedState
            );

        std::cout << "\nMass balance after one step:\n";
        std::cout << reactorModel.massBalanceToString(
            massBalanceAfterStep
        ) << std::endl;
    }

    void runImplicitEulerSingleStepDiagnostic(
        const ReactorState& initializedState,
        const ReducedReactorDae& reducedDae,
        const WellStirredReactorModel& reactorModel
    )
    {
        const ReducedReactorDaeStepResult implicitStep =
            reducedDae.implicitEulerStep(
                initializedState,
                1.0
            );

        printSectionTitle("Reduced DAE implicit Euler diagnostic");
        std::cout << reducedDae.stepResultToString(implicitStep) << std::endl;

        const MassBalanceResult implicitStepMassBalance =
            reactorModel.computeMassBalance(implicitStep.state);

        std::cout << "\nMass balance after implicit Euler step:\n";
        std::cout << reactorModel.massBalanceToString(
            implicitStepMassBalance
        ) << std::endl;
    }

    void runRadauSingleStepDiagnostic(
        const ReactorState& initializedState,
        const ReducedReactorDae& reducedDae,
        const WellStirredReactorModel& reactorModel
    )
    {
        const ReducedReactorDaeStepResult radauStep =
            reducedDae.radauIIA2Step(
                initializedState,
                1.0
            );

        printSectionTitle("Reduced DAE Radau IIA 2-stage diagnostic");
        std::cout << reducedDae.stepResultToString(radauStep) << std::endl;

        const MassBalanceResult radauStepMassBalance =
            reactorModel.computeMassBalance(radauStep.state);

        std::cout << "\nMass balance after Radau IIA 2-stage step:\n";
        std::cout << reactorModel.massBalanceToString(
            radauStepMassBalance
        ) << std::endl;
    }

    void runEnthalpyDiagnostic(
        const ReactorState& initializedState,
        const ThermoPackage& thermo,
        const WellStirredReactorModel& reactorModel
        )
    {
        const WellStirredReactorParameters& parameters =
            reactorModel.parameters();

        const WellStirredReactorEvaluation evaluation =
            reactorModel.evaluate(initializedState);

        const MassBalanceResult massBalance =
            reactorModel.computeMassBalance(initializedState);

        const double reactorTemperatureK =
            initializedState.temperatureC() + 273.15;

        const double inletTemperatureK =
            parameters.inletTemperatureC + 273.15;

        const double inletPressurePa =
            parameters.inletPressureBar * 1.0e5;

        const FlashResult inletFlash =
            thermo.flash(
                inletPressurePa,
                inletTemperatureK,
                parameters.inletComposition
                );

        const double inletMolarEnthalpy =
            inletFlash.molarEnthalpyMixtureJPerKmol;

        const double reactorMolarEnthalpy =
            evaluation.flash.molarEnthalpyMixtureJPerKmol;

        const double outletMolarEnthalpy =
            evaluation.flash.molarEnthalpyVaporJPerKmol;

        const double inletEnthalpyFlow =
            thermo.enthalpy().enthalpyFlowW(
                parameters.inletFlowKmolPerS,
                inletTemperatureK,
                parameters.inletComposition
                );

        const double outletEnthalpyFlow =
            thermo.enthalpy().enthalpyFlowW(
                massBalance.outletFlowKmolPerS,
                reactorTemperatureK,
                evaluation.flash.yVapor
                );

        const double estimatedInventoryEnergy =
            evaluation.totalMolesKmol
            * reactorMolarEnthalpy;

        printSectionTitle("Enthalpy diagnostic");

        std::cout << std::fixed << std::setprecision(6);

        std::cout << "Molar enthalpies:\n";
        std::cout << "  h_inlet   = "
                  << inletMolarEnthalpy
                  << " J/kmol\n";
        std::cout << "  h_reactor = "
                  << reactorMolarEnthalpy
                  << " J/kmol\n";
        std::cout << "  h_outlet  = "
                  << outletMolarEnthalpy
                  << " J/kmol\n";

        std::cout << "\nEnthalpy flows:\n";
        std::cout << "  Hdot_in   = "
                  << inletEnthalpyFlow
                  << " W\n";
        std::cout << "  Hdot_out  = "
                  << outletEnthalpyFlow
                  << " W\n";
        std::cout << "  net Hdot  = "
                  << inletEnthalpyFlow - outletEnthalpyFlow
                  << " W\n";

        std::cout << "\nInventory energy estimate:\n";
        std::cout << "  E_est     = "
                  << estimatedInventoryEnergy
                  << " J\n";
        std::cout << "  E_state   = "
                  << initializedState.energyJ()
                  << " J\n";
        std::cout << "  E_error   = "
                  << initializedState.energyJ() - estimatedInventoryEnergy
                  << " J\n";
    }

    void runEnergyDaeResidualDiagnostic(
        const ReactorState& initializedState,
        const EnergyReactorDae& energyDae,
        const WellStirredReactorModel& reactorModel
        )
    {
        const MassBalanceResult massBalance =
            reactorModel.computeMassBalance(
                initializedState
                );

        const EnergyBalanceResult energyBalance =
            reactorModel.computeEnergyBalance(
                initializedState
                );

        EnergyReactorDaeVariables variables =
            energyDae.variablesFromState(
                initializedState
                );

        EnergyReactorDaeDerivatives derivatives;

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            derivatives.massDerivativesKgPerS[i] =
                massBalance.massDerivativesKgPerS[i];
        }

        derivatives.energyDerivativeW =
            energyBalance.energyDerivativeW;

        const EnergyReactorDaeResidual residual =
            energyDae.computeResidual(
                variables,
                derivatives
                );

        printSectionTitle("Energy DAE residual diagnostic");

        std::cout << energyDae.residualToString(
            residual
            ) << std::endl;

        const numerics::Vector variablesVector =
            energyDae.variablesToVector(
                variables
                );

        const numerics::Vector derivativesVector =
            energyDae.derivativesToVector(
                derivatives
                );

        const numerics::Vector residualVector =
            energyDae.computeResidualVector(
                variablesVector,
                derivativesVector
                );

        printSectionTitle("Energy DAE vector residual diagnostic");

        std::cout << std::fixed << std::setprecision(10);

        for (std::size_t i = 0; i < residualVector.size(); ++i)
        {
            std::cout << "R[" << i << "] = "
                      << residualVector[i]
                      << "\n";
        }
    }

    void runEnergyRadauSingleStepDiagnostic(
        const ReactorState& initializedState,
        const EnergyReactorDae& energyDae,
        const WellStirredReactorModel& reactorModel
        )
    {
        const EnergyReactorDaeStepResult stepResult =
            energyDae.radauIIA3Step(
                initializedState,
                1.0
                );

        printSectionTitle("Energy DAE Radau IIA 3-stage diagnostic");

        std::cout << energyDae.stepResultToString(
            stepResult
            ) << std::endl;

        const MassBalanceResult massBalance =
            reactorModel.computeMassBalance(
                stepResult.state
                );

        const EnergyBalanceResult energyBalance =
            reactorModel.computeEnergyBalance(
                stepResult.state
                );

        std::cout << "\nMass balance after energy Radau step:\n";
        std::cout << reactorModel.massBalanceToString(
            massBalance
            ) << std::endl;

        std::cout << "\nEnergy balance after energy Radau step:\n";
        std::cout << reactorModel.energyBalanceToString(
            energyBalance
            ) << std::endl;
    }

    double mixtureHeatCapacityJPerKmolK(
        const ThermoPackage& thermo,
        const Composition& composition
        )
    {
        const ComponentEnthalpyDataList& data =
            thermo.enthalpy().data();

        double compositionSum = 0.0;
        double cpSum = 0.0;

        for (std::size_t i = 0; i < ComponentCount; ++i)
        {
            compositionSum +=
                composition[i];

            cpSum +=
                composition[i]
                * data[i].cpA;
        }

        if (compositionSum <= 0.0)
        {
            throw std::runtime_error("Cannot calculate mixture Cp for empty composition");
        }

        return cpSum / compositionSum;
    }

    void runClosedHeatingValidation(
        const ReactorState& initializedState,
        const ThermoPackage& thermo,
        const WellStirredReactorParameters& baseParameters
        )
    {
        printSectionTitle("Closed heating validation");

        ReactionModel validationReactionModel(
            thermo.materials()
            );

        validationReactionModel.setForwardRateConstant(0.0);
        validationReactionModel.setReverseRateConstant(0.0);

        WellStirredReactorParameters validationParameters =
            baseParameters;

        validationParameters.inletFlowKmolPerS = 0.0;
        validationParameters.outletValveCoefficientKmolPerSBar = 0.0;
        validationParameters.allowReverseOutletFlow = false;

        validationParameters.heatTransferAreaM2 = 1.0;
        validationParameters.heatTransferCoefficientWPerM2K = 500.0;
        validationParameters.jacketTemperatureC = 110.0;

        WellStirredReactorModel validationModel(
            validationParameters,
            thermo,
            validationReactionModel
            );

        ReactorState validationInitialState =
            validationModel.initializeEnergyFromEnthalpy(
                initializedState
                );

        EnergyReactorDae validationDae(
            validationModel
            );

        EnergyReactorSimulationRunner validationRunner(
            validationDae,
            validationModel
            );

        EnergyReactorSimulationOptions options;
        options.timeStepS = 10.0;
        options.stepCount = 60; // 600 s

        const EnergyReactorSimulationResult result =
            validationRunner.runRadauIIA3(
                validationInitialState,
                options
                );

        if (!result.converged)
        {
            std::cout << "Validation simulation failed: "
                      << result.message
                      << "\n";
            return;
        }

        const ReactorState& finalState =
            result.finalState();

        const WellStirredReactorEvaluation initialEvaluation =
            validationModel.evaluate(
                validationInitialState
                );

        const Composition initialComposition =
            validationInitialState.moleFractions();

        const double totalMolesKmol =
            initialEvaluation.totalMolesKmol;

        const double cpMix =
            mixtureHeatCapacityJPerKmolK(
                thermo,
                initialComposition
                );

        const double heatConductanceWPerK =
            validationParameters.heatTransferCoefficientWPerM2K
            * validationParameters.heatTransferAreaM2;

        const double initialTemperatureC =
            validationInitialState.temperatureC();

        const double jacketTemperatureC =
            validationParameters.jacketTemperatureC;

        const double finalTimeS =
            options.timeStepS * options.stepCount;

        const double timeConstantS =
            totalMolesKmol
            * cpMix
            / heatConductanceWPerK;

        const double analyticalTemperatureC =
            jacketTemperatureC
            - (
                  jacketTemperatureC
                  - initialTemperatureC
                  )
                  * std::exp(
                      -finalTimeS / timeConstantS
                      );

        const double numericalTemperatureC =
            finalState.temperatureC();

        double maxInventoryEnergyResidualJ = 0.0;
        double maxNewtonResidualNorm = 0.0;

        for (const EnergyReactorSimulationPoint& point : result.points)
        {
            maxInventoryEnergyResidualJ =
                std::max(
                    maxInventoryEnergyResidualJ,
                    std::abs(point.energyBalance.energyResidualJ)
                    );

            maxNewtonResidualNorm =
                std::max(
                    maxNewtonResidualNorm,
                    point.residualNorm
                    );
        }

        std::cout << std::fixed << std::setprecision(10);

        std::cout << "Validation case:\n";
        std::cout << "  inlet flow       = 0 kmol/s\n";
        std::cout << "  outlet valve     = closed\n";
        std::cout << "  reaction rates   = 0\n";
        std::cout << "  U                = "
                  << validationParameters.heatTransferCoefficientWPerM2K
                  << " W/(m2*K)\n";
        std::cout << "  A                = "
                  << validationParameters.heatTransferAreaM2
                  << " m2\n";
        std::cout << "  T_jacket         = "
                  << jacketTemperatureC
                  << " C\n";

        std::cout << "\nAnalytical check:\n";
        std::cout << "  total moles N    = "
                  << totalMolesKmol
                  << " kmol\n";
        std::cout << "  Cp_mix           = "
                  << cpMix
                  << " J/(kmol*K)\n";
        std::cout << "  tau              = "
                  << timeConstantS
                  << " s\n";
        std::cout << "  final time       = "
                  << finalTimeS
                  << " s\n";

        std::cout << "\nTemperature comparison:\n";
        std::cout << "  T_initial        = "
                  << initialTemperatureC
                  << " C\n";
        std::cout << "  T_numerical      = "
                  << numericalTemperatureC
                  << " C\n";
        std::cout << "  T_analytical     = "
                  << analyticalTemperatureC
                  << " C\n";
        std::cout << "  absolute error   = "
                  << std::abs(numericalTemperatureC - analyticalTemperatureC)
                  << " C\n";

        std::cout << "\nDAE consistency:\n";
        std::cout << "  max |E - N*h|    = "
                  << maxInventoryEnergyResidualJ
                  << " J\n";
        std::cout << "  max residual norm = "
                  << maxNewtonResidualNorm
                  << "\n";

        std::cout << "\nFinal validation state:\n";
        std::cout << finalState.toString() << "\n";
    }

}

int main(int argc, char* argv[])
{
    ReactorState initialState;
    ThermoPackage thermo;

    ReactionModel reactionModel(thermo.materials());
    reactionModel.setForwardRateConstant(100.0);
    reactionModel.setReverseRateConstant(100.0);

    WellStirredReactorParameters reactorParameters;
    reactorParameters.volumeM3 = 1.0;
    reactorParameters.allowReverseOutletFlow = false;

    WellStirredReactorModel reactorModel(
        reactorParameters,
        thermo,
        reactionModel
    );

    ReducedReactorDae reducedDae(reactorModel);

    EnergyReactorDae energyDae(
        reactorModel
        );

    if (PrintInitialState)
    {
        printSectionTitle("Initial reactor state");
        std::cout << initialState.toString() << std::endl;
    }

    if (PrintMaterials)
    {
        printMaterials(thermo);
    }

    if (RunFlashDiagnostics)
    {
        runFlashDiagnostics(initialState, thermo);
    }

    if (RunReactorModelDiagnostic)
    {
        const WellStirredReactorEvaluation reactorEvaluation =
            reactorModel.evaluate(initialState);

        printSectionTitle("Reactor model diagnostic");
        std::cout << reactorModel.evaluationToString(
            reactorEvaluation
        ) << std::endl;
    }

    const PressureInitializationResult pressureInitialization =
        reactorModel.initializePressureForVolumeNewton(
            initialState,
            InitialPressureSearchMinPa,
            InitialPressureSearchMaxPa,
            InitialVolumeToleranceM3,
            InitialPressureNewtonMaxIterations,
            true
        );

    if (RunPressureInitializationDiagnostic)
    {
        printPressureInitializationSummary(pressureInitialization);
    }

    if (!pressureInitialization.converged)
    {
        std::cerr << "Pressure initialization failed. Stop simulation.\n";
        return 1;
    }

    const ReactorState initializedState =
        reactorModel.initializeEnergyFromEnthalpy(
            pressureInitialization.initializedState
            );

    if (RunMassBalanceDiagnostic)
    {
        const MassBalanceResult initialMassBalance =
            reactorModel.computeMassBalance(initializedState);

        printSectionTitle("Mass balance diagnostic");
        std::cout << reactorModel.massBalanceToString(
            initialMassBalance
            ) << std::endl;
    }

    if (RunEnthalpyDiagnostic)
    {
        runEnthalpyDiagnostic(
            initializedState,
            thermo,
            reactorModel
            );
    }

    if (RunEnergyBalanceDiagnostic)
    {
        const EnergyBalanceResult energyBalance =
            reactorModel.computeEnergyBalance(initializedState);

        printSectionTitle("Energy balance diagnostic");

        std::cout << reactorModel.energyBalanceToString(
            energyBalance
            ) << std::endl;
    }

    if (RunEnergyDaeResidualDiagnostic)
    {
        runEnergyDaeResidualDiagnostic(
            initializedState,
            energyDae,
            reactorModel
            );
    }

    if (RunEnergyRadauSingleStepDiagnostic)
    {
        runEnergyRadauSingleStepDiagnostic(
            initializedState,
            energyDae,
            reactorModel
            );
    }

    if (RunDaeResidualDiagnostic)
    {
        const MassBalanceResult initialMassBalance =
            reactorModel.computeMassBalance(initializedState);

        runDaeResidualDiagnostics(
            initializedState,
            reducedDae,
            initialMassBalance
            );
    }

    if (RunNewtonSimpleDiagnostic)
    {
        runNewtonSimpleDiagnostic();
    }

    if (RunExplicitEulerSingleStepDiagnostic)
    {
        runExplicitEulerSingleStepDiagnostic(
            initializedState,
            reactorModel
        );
    }

    if (RunExplicitEulerSimulationDiagnostic)
    {
        (void)runExplicitEulerMassSimulation(
            initializedState,
            reactorModel,
            0.005,
            200,
            20
        );
    }

    if (RunImplicitEulerSingleStepDiagnostic)
    {
        runImplicitEulerSingleStepDiagnostic(
            initializedState,
            reducedDae,
            reactorModel
        );
    }

    if (RunImplicitEulerSimulationDiagnostic)
    {
        (void)runImplicitEulerSimulation(
            initializedState,
            reducedDae,
            reactorModel,
            MainTimeStepS,
            MainStepCount,
            MainPrintEvery
        );
    }

    if (RunRadauSingleStepDiagnostic)
    {
        runRadauSingleStepDiagnostic(
            initializedState,
            reducedDae,
            reactorModel
        );
    }

    ReactorState finalState = initializedState;

    ReactorSimulationResult simulationResult;

    if (RunRadauSimulationDiagnostic)
    {
        ReactorSimulationRunner simulationRunner(
            reducedDae,
            reactorModel
            );

        ReactorSimulationOptions simulationOptions;
        simulationOptions.timeStepS = MainTimeStepS;
        simulationOptions.stepCount = MainStepCount;

        simulationResult =
            simulationRunner.runRadauIIA3(
                initializedState,
                simulationOptions
                );

        printSectionTitle("Radau IIA 3-stage DAE simulation diagnostic");

        std::cout << "dt = "
                  << simulationOptions.timeStepS
                  << " s\n";

        std::cout << "steps = "
                  << simulationOptions.stepCount
                  << "\n\n";

        std::cout << simulationRunner.resultTableToString(
            simulationResult,
            MainPrintEvery
            ) << std::endl;

        if (!simulationResult.converged)
        {
            std::cout << "Simulation failed: "
                      << simulationResult.message
                      << std::endl;
        }

        if (!simulationResult.points.empty())
        {
            finalState =
                simulationResult.finalState();
        }
    }

    if (RunEnergyRadauSimulationDiagnostic)
    {
        EnergyReactorSimulationRunner energySimulationRunner(
            energyDae,
            reactorModel
            );

        EnergyReactorSimulationOptions energySimulationOptions;
        energySimulationOptions.timeStepS = MainTimeStepS;
        energySimulationOptions.stepCount = MainStepCount;

        EnergyReactorSimulationResult energySimulationResult =
            energySimulationRunner.runRadauIIA3(
                initializedState,
                energySimulationOptions
                );

        printSectionTitle(
            "Energy Radau IIA 3-stage DAE simulation diagnostic"
            );

        std::cout << "dt = "
                  << energySimulationOptions.timeStepS
                  << " s\n";

        std::cout << "steps = "
                  << energySimulationOptions.stepCount
                  << "\n\n";

        std::cout << energySimulationRunner.resultTableToString(
            energySimulationResult,
            MainPrintEvery
            ) << std::endl;

        if (!energySimulationResult.converged)
        {
            std::cout << "Energy simulation failed: "
                      << energySimulationResult.message
                      << std::endl;
        }

        if (!energySimulationResult.points.empty())
        {
            finalState =
                energySimulationResult.finalState();
        }
    }

    if (RunClosedHeatingValidation)
    {
        runClosedHeatingValidation(
            initializedState,
            thermo,
            reactorParameters
            );
    }

    if (PrintFinalState)
    {
        printSectionTitle("Final reactor state");
        std::cout << finalState.toString() << std::endl;
    }

    QApplication app(argc, argv);

    MainWindow window;
    window.show();

    return app.exec();
}
