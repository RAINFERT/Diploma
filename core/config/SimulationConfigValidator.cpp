#include "SimulationConfigValidator.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace {

void requireFinite(
    const double value,
    const char* name
    )
{
    if (!std::isfinite(value)) {
        std::ostringstream message;
        message << "Config value must be finite: " << name;
        throw std::runtime_error(message.str());
    }
}

void requireRange(
    const double value,
    const double minValue,
    const double maxValue,
    const char* name
    )
{
    requireFinite(value, name);

    if (value < minValue || value > maxValue) {
        std::ostringstream message;
        message
            << "Config value out of range: " << name
            << ", value = " << value
            << ", allowed range = [" << minValue
            << ", " << maxValue << "]";

        throw std::runtime_error(message.str());
    }
}

void requireIntRange(
    const int value,
    const int minValue,
    const int maxValue,
    const char* name
    )
{
    if (value < minValue || value > maxValue) {
        std::ostringstream message;
        message
            << "Config integer value out of range: " << name
            << ", value = " << value
            << ", allowed range = [" << minValue
            << ", " << maxValue << "]";

        throw std::runtime_error(message.str());
    }
}

void requireMinLessThanMax(
    const double minValue,
    const double maxValue,
    const char* minName,
    const char* maxName
    )
{
    requireFinite(minValue, minName);
    requireFinite(maxValue, maxName);

    if (minValue >= maxValue) {
        std::ostringstream message;
        message
            << "Invalid limits: " << minName
            << " must be less than " << maxName
            << ". Got " << minValue
            << " >= " << maxValue;

        throw std::runtime_error(message.str());
    }
}

void validateLimitsBlock(const LimitsConfig& limits)
{
    requireMinLessThanMax(
        limits.minComponentMassKg,
        limits.maxComponentMassKg,
        "limits.minComponentMassKg",
        "limits.maxComponentMassKg"
        );

    requireMinLessThanMax(
        limits.minTemperatureC,
        limits.maxTemperatureC,
        "limits.minTemperatureC",
        "limits.maxTemperatureC"
        );

    requireMinLessThanMax(
        limits.minPressureBar,
        limits.maxPressureBar,
        "limits.minPressureBar",
        "limits.maxPressureBar"
        );

    requireMinLessThanMax(
        limits.minFlowKmolPerS,
        limits.maxFlowKmolPerS,
        "limits.minFlowKmolPerS",
        "limits.maxFlowKmolPerS"
        );

    requireMinLessThanMax(
        limits.minReactorVolumeM3,
        limits.maxReactorVolumeM3,
        "limits.minReactorVolumeM3",
        "limits.maxReactorVolumeM3"
        );

    requireMinLessThanMax(
        limits.minTimeStepS,
        limits.maxTimeStepS,
        "limits.minTimeStepS",
        "limits.maxTimeStepS"
        );

    if (limits.minStepCount >= limits.maxStepCount) {
        throw std::runtime_error(
            "Invalid limits: limits.minStepCount must be less than limits.maxStepCount"
            );
    }

    requireRange(
        limits.compositionSumTolerance,
        0.0,
        1.0e-2,
        "limits.compositionSumTolerance"
        );
}

void requireNonEmptyString(
    const std::string& value,
    const char* name
    )
{
    if (value.empty()) {
        std::ostringstream message;
        message << "Config string must not be empty: " << name;
        throw std::runtime_error(message.str());
    }
}

void validateOutputBlock(const OutputConfig& output)
{
    if (!output.enabled) {
        return;
    }

    requireNonEmptyString(
        output.directory,
        "output.directory"
        );

    if (output.writeResultTable) {
        requireNonEmptyString(
            output.resultTableFile,
            "output.resultTableFile"
            );
    }

    if (output.writeInitialState) {
        requireNonEmptyString(
            output.initialStateFile,
            "output.initialStateFile"
            );
    }

    if (output.writePressureInitialization) {
        requireNonEmptyString(
            output.pressureInitializationFile,
            "output.pressureInitializationFile"
            );
    }

    if (output.writeFinalState) {
        requireNonEmptyString(
            output.finalStateFile,
            "output.finalStateFile"
            );
    }
}

void validateComposition(
    const Composition& composition,
    const LimitsConfig& limits,
    const char* name
    )
{
    double sum = 0.0;

    for (std::size_t i = 0; i < ComponentCount; ++i) {
        const double value = composition[i];

        requireFinite(value, name);

        if (value < 0.0) {
            throw std::runtime_error("Composition contains negative component fraction");
        }

        if (value > 1.0) {
            throw std::runtime_error("Composition contains component fraction greater than 1");
        }

        sum += value;
    }

    if (std::abs(sum - 1.0) > limits.compositionSumTolerance) {
        std::ostringstream message;
        message
            << "Composition sum is not 1.0 for " << name
            << ". Actual sum = " << sum
            << ", tolerance = " << limits.compositionSumTolerance;

        throw std::runtime_error(message.str());
    }
}

} // namespace

void SimulationConfigValidator::validate(const SimulationConfig& config)
{
    const LimitsConfig& limits = config.limits;

    validateLimitsBlock(limits);
    validateOutputBlock(config.output);

    if (config.schemaVersion != 1) {
        throw std::runtime_error("Unsupported config schemaVersion");
    }

    if (config.caseName.empty()) {
        throw std::runtime_error("Config caseName must not be empty");
    }

    requireNonEmptyString(
        config.dataSources.componentPropertiesCsv,
        "dataSources.componentPropertiesCsv"
        );

    if (config.components.size() != ComponentCount) {
        throw std::runtime_error(
            "At stage 1 config must contain exactly 3 components"
            );
    }

    if (config.componentKeys.size() != config.components.size()) {
        throw std::runtime_error(
            "Internal config error: componentKeys size differs from components size"
            );
    }

    requireRange(
        config.feed.pressureBar,
        limits.minPressureBar,
        limits.maxPressureBar,
        "feed.pressureBar"
        );

    requireRange(
        config.feed.temperatureC,
        limits.minTemperatureC,
        limits.maxTemperatureC,
        "feed.temperatureC"
        );

    requireRange(
        config.feed.flowKmolPerS,
        limits.minFlowKmolPerS,
        limits.maxFlowKmolPerS,
        "feed.flowKmolPerS"
        );

    validateComposition(
        config.feed.composition,
        limits,
        "feed.composition"
        );

    requireRange(
        config.reactor.volumeM3,
        limits.minReactorVolumeM3,
        limits.maxReactorVolumeM3,
        "reactor.volumeM3"
        );

    requireRange(
        config.reactor.outlet.valveCoefficientKmolPerSBar,
        0.0,
        limits.maxFlowKmolPerS,
        "reactor.outlet.valveCoefficientKmolPerSBar"
        );

    requireRange(
        config.reactor.outlet.referencePressureBar,
        limits.minPressureBar,
        limits.maxPressureBar,
        "reactor.outlet.referencePressureBar"
        );

    requireRange(
        config.reactor.heatTransfer.areaM2,
        0.0,
        1.0e9,
        "reactor.heatTransfer.areaM2"
        );

    requireRange(
        config.reactor.heatTransfer.coefficientWPerM2K,
        0.0,
        1.0e9,
        "reactor.heatTransfer.coefficientWPerM2K"
        );

    requireRange(
        config.reactor.heatTransfer.jacketTemperatureC,
        limits.minTemperatureC,
        limits.maxTemperatureC,
        "reactor.heatTransfer.jacketTemperatureC"
        );

    requireRange(
        config.initialState.temperatureC,
        limits.minTemperatureC,
        limits.maxTemperatureC,
        "initialState.temperatureC"
        );

    requireRange(
        config.initialState.pressureBar,
        limits.minPressureBar,
        limits.maxPressureBar,
        "initialState.pressureBar"
        );

    for (std::size_t i = 0; i < ComponentCount; ++i) {
        requireRange(
            config.initialState.componentMassesKg[i],
            limits.minComponentMassKg,
            limits.maxComponentMassKg,
            "initialState.componentMassesKg"
            );
    }

    requireMinLessThanMax(
        config.initialState.pressureInitialization.pressureMinPa,
        config.initialState.pressureInitialization.pressureMaxPa,
        "initialState.pressureInitialization.pressureMinPa",
        "initialState.pressureInitialization.pressureMaxPa"
        );

    requireRange(
        config.initialState.pressureInitialization.volumeToleranceM3,
        1.0e-15,
        1.0,
        "initialState.pressureInitialization.volumeToleranceM3"
        );

    requireIntRange(
        config.initialState.pressureInitialization.newtonMaxIterations,
        1,
        10000,
        "initialState.pressureInitialization.newtonMaxIterations"
        );

    requireRange(
        config.reactions.forwardRateConstant1PerS,
        0.0,
        1.0e12,
        "reactions.forwardRateConstant1PerS"
        );

    requireRange(
        config.reactions.reverseRateConstant1PerS,
        0.0,
        1.0e12,
        "reactions.reverseRateConstant1PerS"
        );

    requireRange(
        config.solver.timeStepS,
        limits.minTimeStepS,
        limits.maxTimeStepS,
        "solver.timeStepS"
        );

    requireIntRange(
        config.solver.stepCount,
        limits.minStepCount,
        limits.maxStepCount,
        "solver.stepCount"
        );

    requireIntRange(
        config.solver.printEvery,
        1,
        config.solver.stepCount,
        "solver.printEvery"
        );
}
