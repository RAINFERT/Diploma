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
    if (printEvery <= 0) {
        throw std::invalid_argument(
            "Print interval must be positive"
            );
    }

    if (result.points.empty()) {
        return {};
    }

    const ReactorState& firstState =
        result.points.front().state;

    const std::size_t componentCount =
        firstState.componentCount();

    if (componentCount == 0) {
        throw std::runtime_error(
            "Cannot print result table: reactor state has no components"
            );
    }

    std::ostringstream out;
    out << std::fixed << std::setprecision(6);

    out
        << "time_s"
        << "\tP_bar"
        << "\tT_C"
        << "\tbeta"
        << "\tFout";

    for (std::size_t i = 0; i < componentCount; ++i) {
        out
            << "\tM_"
            << firstState.componentName(i)
            << "_kg";
    }

    out << "\tE_J";

    for (std::size_t i = 0; i < componentCount; ++i) {
        out
            << "\tdM_"
            << firstState.componentName(i)
            << "_dt";
    }

    out
        << "\tdE_dt_W"
        << "\tQdot_W"
        << "\tHdot_in_W"
        << "\tHdot_out_W"
        << "\tNewton_it"
        << "\n";

    for (std::size_t index = 0; index < result.points.size(); ++index) {
        const bool shouldPrint =
            index % static_cast<std::size_t>(printEvery) == 0 ||
            index + 1 == result.points.size();

        if (!shouldPrint) {
            continue;
        }

        const EnergyReactorSimulationPoint& point =
            result.points[index];

        const ReactorState& state =
            point.state;

        if (state.componentCount() != componentCount) {
            throw std::runtime_error(
                "Cannot print result table: component count changed during simulation"
                );
        }

        const WellStirredReactorEvaluation& evaluation =
            point.evaluation;

        const MassBalanceResult& massBalance =
            point.massBalance;

        const EnergyBalanceResult& energyBalance =
            point.energyBalance;

        if (massBalance.massDerivativesKgPerS.size() != componentCount) {
            throw std::runtime_error(
                "Cannot print result table: mass derivative vector size differs from component count"
                );
        }

        out
            << point.timeS
            << "\t" << state.pressureBar()
            << "\t" << state.temperatureC()
            << "\t" << evaluation.flash.beta
            << "\t" << massBalance.outletFlowKmolPerS;

        for (std::size_t i = 0; i < componentCount; ++i) {
            out
                << "\t"
                << state.massKg(i);
        }

        out
            << "\t"
            << state.energyJ();

        for (std::size_t i = 0; i < componentCount; ++i) {
            out
                << "\t"
                << massBalance.massDerivativesKgPerS[i];
        }

        out
            << "\t" << energyBalance.energyDerivativeW
            << "\t" << energyBalance.heatTransferW
            << "\t" << energyBalance.inletEnthalpyFlowW
            << "\t" << energyBalance.outletEnthalpyFlowW
            << "\t" << point.newtonIterations
            << "\n";
    }

    return out.str();
}
