#include "DynamicEnergyReactorSimulation.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>

DynamicEnergyReactorSimulationRunner::DynamicEnergyReactorSimulationRunner(
    const DynamicEnergyReactorDae& dae,
    const DynamicWellStirredReactorModel& reactorModel
    )
    : dae_(dae),
    reactorModel_(reactorModel)
{
}

DynamicEnergyReactorSimulationResult
DynamicEnergyReactorSimulationRunner::runRadauIIA3(
    const DynamicVector& initialVariables,
    const DynamicEnergyReactorSimulationOptions& options
    ) const {
    if (options.timeStepS <= 0.0) {
        throw std::runtime_error(
            "DynamicEnergyReactorSimulationRunner: time step must be positive"
            );
    }

    if (options.stepCount <= 0) {
        throw std::runtime_error(
            "DynamicEnergyReactorSimulationRunner: step count must be positive"
            );
    }

    DynamicEnergyReactorSimulationResult result;

    DynamicEnergyReactorRadau radau(
        dae_,
        reactorModel_
        );

    DynamicVector variables =
        initialVariables;

    for (int step = 0; step <= options.stepCount; ++step) {
        DynamicEnergyReactorDaeEvaluation evaluation =
            dae_.evaluate(variables);

        DynamicMassBalanceResult massBalance =
            reactorModel_.computeMassBalance(
                evaluation.pressurePa,
                evaluation.massesKg,
                evaluation.gibbsOde.yVapor
                );

        DynamicEnergyBalanceResult energyBalance =
            reactorModel_.computeEnergyBalance(
                evaluation.pressurePa,
                evaluation.temperatureK,
                massBalance.outletMolarFlowKmolPerS,
                evaluation.gibbsOde.yVapor
                );

        DynamicEnergyReactorSimulationPoint point;
        point.timeS =
            step * options.timeStepS;

        point.variables =
            variables;

        point.evaluation =
            evaluation;

        point.massBalance =
            massBalance;

        point.energyBalance =
            energyBalance;

        result.points.push_back(point);

        if (step == options.stepCount) {
            break;
        }

        DynamicRadauStepResult stepResult =
            radau.radauIIA3Step(
                variables,
                options.timeStepS,
                options.radauOptions
                );

        if (!stepResult.converged) {
            result.converged = false;
            result.message =
                "Dynamic Radau failed at step " +
                std::to_string(step + 1) +
                ": " +
                stepResult.message;

            DynamicEnergyReactorSimulationPoint failedPoint;
            failedPoint.timeS =
                (step + 1) * options.timeStepS;
            failedPoint.variables =
                stepResult.newVariables;
            failedPoint.converged =
                false;
            failedPoint.newtonIterations =
                stepResult.iterations;
            failedPoint.residualNorm =
                stepResult.residualNorm;

            result.points.push_back(failedPoint);
            return result;
        }

        variables =
            stepResult.newVariables;

        result.points.back().newtonIterations =
            stepResult.iterations;

        result.points.back().residualNorm =
            stepResult.residualNorm;
    }

    result.converged = true;
    result.message = "Dynamic Radau simulation converged";

    return result;
}

std::string DynamicEnergyReactorSimulationRunner::resultTableToString(
    const DynamicEnergyReactorSimulationResult& result,
    const ComponentSet& components,
    int printEvery
    ) const {
    if (printEvery <= 0) {
        throw std::runtime_error(
            "DynamicEnergyReactorSimulationRunner: printEvery must be positive"
            );
    }

    std::ostringstream stream;

    stream << std::fixed << std::setprecision(6);

    stream
        << "time_s"
        << "\tP_bar"
        << "\tT_K"
        << "\tbetaV"
        << "\tFout_kmol_s"
        << "\tdE_dt_W"
        << "\tQdot_W"
        << "\tHdot_in_W"
        << "\tHdot_out_W"
        << "\tE_inventory_J"
        << "\thMix_J_kmol"
        << "\thL_total_J_kmol"
        << "\thV_total_J_kmol"
        << "\thL_dep_J_kmol"
        << "\thV_dep_J_kmol";

    for (const ComponentProperties& component : components.components) {
        stream << "\tM_" << component.name << "_kg";
    }

    stream
        << "\tvolumeResidual"
        << "\tflashResidual"
        << "\tNewton_it"
        << "\tNewton_residual"
        << "\n";

    for (std::size_t pointIndex = 0;
         pointIndex < result.points.size();
         ++pointIndex) {

        if (static_cast<int>(pointIndex) % printEvery != 0 &&
            pointIndex + 1 != result.points.size()) {
            continue;
        }

        const DynamicEnergyReactorSimulationPoint& point =
            result.points[pointIndex];

        if (!point.converged) {
            stream
                << point.timeS
                << "\tFAILED\n";
            continue;
        }

        const DynamicEnergyReactorDaeEvaluation& evaluation =
            point.evaluation;

        stream
            << point.timeS
            << "\t" << evaluation.pressurePa / 1.0e5
            << "\t" << evaluation.temperatureK
            << "\t" << evaluation.gibbsOde.betaVapor
            << "\t" << point.massBalance.outletMolarFlowKmolPerS
            << "\t" << point.energyBalance.energyRateW
            << "\t" << point.energyBalance.heatTransferW
            << "\t" << point.energyBalance.inletEnthalpyFlowW
            << "\t" << point.energyBalance.outletEnthalpyFlowW
            << "\t" << evaluation.calculatedInventoryEnergyJ
            << "\t" << evaluation.twoPhaseMolarEnthalpyJPerKmol
            << "\t" << evaluation.liquidTotalEnthalpyJPerKmol
            << "\t" << evaluation.vaporTotalEnthalpyJPerKmol
            << "\t" << evaluation.liquidDepartureEnthalpyJPerKmol
            << "\t" << evaluation.vaporDepartureEnthalpyJPerKmol;

        for (double massKg : evaluation.massesKg) {
            stream << "\t" << massKg;
        }

        stream
            << "\t" << evaluation.volumeResidual
            << "\t" << evaluation.gibbsOde.maxFugacityResidual
            << "\t" << point.newtonIterations
            << "\t" << point.residualNorm
            << "\n";
    }

    return stream.str();
}
