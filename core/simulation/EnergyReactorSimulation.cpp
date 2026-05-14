#include "EnergyReactorSimulation.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>

const ReactorState& EnergyReactorSimulationResult::finalState() const
{
    if (points.empty())
    {
        throw std::runtime_error("Energy simulation result has no points");
    }

    return points.back().state;
}

EnergyReactorSimulationRunner::EnergyReactorSimulationRunner(
    const EnergyReactorDae& dae,
    const WellStirredReactorModel& reactorModel
    )
    : dae_(dae),
    reactorModel_(reactorModel)
{
}

EnergyReactorSimulationPoint EnergyReactorSimulationRunner::makePoint(
    double timeS,
    const ReactorState& state,
    int newtonIterations,
    double residualNorm,
    bool converged,
    const std::string& message
    ) const
{
    EnergyReactorSimulationPoint point;

    point.timeS = timeS;
    point.state = state;

    point.evaluation =
        reactorModel_.evaluate(state);

    point.massBalance =
        reactorModel_.computeMassBalance(state);

    point.energyBalance =
        reactorModel_.computeEnergyBalance(state);

    point.newtonIterations = newtonIterations;
    point.residualNorm = residualNorm;
    point.converged = converged;
    point.message = message;

    return point;
}

EnergyReactorSimulationResult EnergyReactorSimulationRunner::runRadauIIA3(
    const ReactorState& initialState,
    const EnergyReactorSimulationOptions& options
    ) const
{
    if (options.timeStepS <= 0.0)
    {
        throw std::invalid_argument("Energy simulation time step must be positive");
    }

    if (options.stepCount <= 0)
    {
        throw std::invalid_argument("Energy simulation step count must be positive");
    }

    EnergyReactorSimulationResult result;

    ReactorState state =
        initialState;

    result.points.push_back(
        makePoint(
            0.0,
            state,
            0,
            0.0,
            true,
            "initial state"
            )
        );

    for (int step = 0; step < options.stepCount; ++step)
    {
        EnergyReactorDaeStepResult stepResult =
            dae_.radauIIA3Step(
                state,
                options.timeStepS
                );

        const double timeS =
            (step + 1) * options.timeStepS;

        if (!stepResult.converged)
        {
            result.converged = false;
            result.message =
                "Energy Radau IIA 3-stage failed at step "
                + std::to_string(step + 1)
                + ": "
                + stepResult.message;

            result.points.push_back(
                makePoint(
                    timeS,
                    state,
                    stepResult.iterations,
                    stepResult.residualNorm,
                    false,
                    stepResult.message
                    )
                );

            return result;
        }

        state =
            stepResult.state;

        result.points.push_back(
            makePoint(
                timeS,
                state,
                stepResult.iterations,
                stepResult.residualNorm,
                true,
                stepResult.message
                )
            );
    }

    result.converged = true;
    result.message = "Energy simulation completed successfully";

    return result;
}

std::string EnergyReactorSimulationRunner::resultTableToString(
    const EnergyReactorSimulationResult& result,
    int printEvery
    ) const
{
    if (printEvery <= 0)
    {
        throw std::invalid_argument("Print interval must be positive");
    }

    std::ostringstream out;

    out << std::fixed << std::setprecision(6);

    out << "time_s"
        << "\tP_bar"
        << "\tT_C"
        << "\tbeta"
        << "\tFout"
        << "\tM_C2H6_kg"
        << "\tM_C5H12_kg"
        << "\tM_H2O_kg"
        << "\tE_J"
        << "\tdM_C2H6_dt"
        << "\tdM_C5H12_dt"
        << "\tdM_H2O_dt"
        << "\tdE_dt_W"
        << "\tQdot_W"
        << "\tHdot_in_W"
        << "\tHdot_out_W"
        << "\tNewton_it"
        << "\n";

    for (std::size_t index = 0; index < result.points.size(); ++index)
    {
        const bool shouldPrint =
            index % static_cast<std::size_t>(printEvery) == 0
            || index + 1 == result.points.size();

        if (!shouldPrint)
        {
            continue;
        }

        const EnergyReactorSimulationPoint& point =
            result.points[index];

        const ReactorState& state =
            point.state;

        const WellStirredReactorEvaluation& evaluation =
            point.evaluation;

        const MassBalanceResult& massBalance =
            point.massBalance;

        const EnergyBalanceResult& energyBalance =
            point.energyBalance;

        out << point.timeS
            << "\t" << state.pressureBar()
            << "\t" << state.temperatureC()
            << "\t" << evaluation.flash.beta
            << "\t" << massBalance.outletFlowKmolPerS
            << "\t" << state.massKg(Component::C2H6)
            << "\t" << state.massKg(Component::C5H12)
            << "\t" << state.massKg(Component::H2O)
            << "\t" << state.energyJ()
            << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::C2H6)]
            << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::C5H12)]
            << "\t" << massBalance.massDerivativesKgPerS[componentIndex(Component::H2O)]
            << "\t" << energyBalance.energyDerivativeW
            << "\t" << energyBalance.heatTransferW
            << "\t" << energyBalance.inletEnthalpyFlowW
            << "\t" << energyBalance.outletEnthalpyFlowW
            << "\t" << point.newtonIterations
            << "\n";
    }

    return out.str();
}
