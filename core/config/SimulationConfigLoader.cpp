#include "SimulationConfigLoader.h"
#include "SimulationConfigValidator.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QString>

#include <cmath>
#include <stdexcept>

namespace {

std::runtime_error configError(const QString& message)
{
    return std::runtime_error(message.toStdString());
}

QJsonObject requireObject(const QJsonObject& parent, const QString& key)
{
    const QJsonValue value = parent.value(key);

    if (!value.isObject()) {
        throw configError("Expected object: " + key);
    }

    return value.toObject();
}

QJsonArray requireArray(const QJsonObject& parent, const QString& key)
{
    const QJsonValue value = parent.value(key);

    if (!value.isArray()) {
        throw configError("Expected array: " + key);
    }

    return value.toArray();
}

double readDouble(
    const QJsonObject& object,
    const QString& key,
    double defaultValue
    )
{
    const QJsonValue value = object.value(key);

    if (value.isUndefined()) {
        return defaultValue;
    }

    if (!value.isDouble()) {
        throw configError("Expected numeric value: " + key);
    }

    return value.toDouble();
}

int readInt(
    const QJsonObject& object,
    const QString& key,
    int defaultValue
    )
{
    const QJsonValue value = object.value(key);

    if (value.isUndefined()) {
        return defaultValue;
    }

    if (!value.isDouble()) {
        throw configError("Expected integer value: " + key);
    }

    const double asDouble = value.toDouble();
    const int asInt = static_cast<int>(asDouble);

    if (std::abs(asDouble - static_cast<double>(asInt)) > 1.0e-12) {
        throw configError("Expected integer value: " + key);
    }

    return asInt;
}

bool readBool(
    const QJsonObject& object,
    const QString& key,
    bool defaultValue
    )
{
    const QJsonValue value = object.value(key);

    if (value.isUndefined()) {
        return defaultValue;
    }

    if (!value.isBool()) {
        throw configError("Expected boolean value: " + key);
    }

    return value.toBool();
}

std::string readString(
    const QJsonObject& object,
    const QString& key,
    const std::string& defaultValue
    )
{
    const QJsonValue value = object.value(key);

    if (value.isUndefined()) {
        return defaultValue;
    }

    if (!value.isString()) {
        throw configError("Expected string value: " + key);
    }

    return value.toString().toStdString();
}

Component componentFromString(const QString& name)
{
    if (name == "C2H6") {
        return Component::C2H6;
    }

    if (name == "C5H12") {
        return Component::C5H12;
    }

    if (name == "H2O") {
        return Component::H2O;
    }

    throw configError("Unknown component: " + name);
}

ThermoFlashMode flashModeFromString(const QString& name)
{
    if (name == "RachfordRice") {
        return ThermoFlashMode::RachfordRice;
    }

    if (name == "GibbsRadauHybrid") {
        return ThermoFlashMode::GibbsRadauHybrid;
    }

    throw configError("Unknown thermo.flash_mode: " + name);
}

double readArrayNumber(
    const QJsonArray& array,
    const int index,
    const QString& arrayName
    )
{
    const QJsonValue value = array.at(index);

    if (!value.isDouble()) {
        throw configError(
            "Expected numeric value in " +
            arrayName +
            " at index " +
            QString::number(index)
            );
    }

    const double number = value.toDouble();

    if (!std::isfinite(number)) {
        throw configError(
            "Expected finite numeric value in " +
            arrayName +
            " at index " +
            QString::number(index)
            );
    }

    return number;
}

void requireArraySize(
    const QJsonArray& array,
    const std::size_t expectedSize,
    const QString& arrayName
    )
{
    if (static_cast<std::size_t>(array.size()) != expectedSize) {
        throw configError(
            "Invalid array size for " +
            arrayName +
            ". Expected " +
            QString::number(static_cast<int>(expectedSize)) +
            ", got " +
            QString::number(array.size())
            );
    }
}

Composition readCompositionArray(
    const QJsonArray& array,
    const std::vector<Component>& components
    )
{
    requireArraySize(
        array,
        components.size(),
        "feed.composition"
        );

    Composition composition{};

    for (std::size_t i = 0; i < components.size(); ++i) {
        const Component component = components[i];

        const double value = readArrayNumber(
            array,
            static_cast<int>(i),
            "feed.composition"
            );

        if (value < 0.0) {
            throw configError(
                "Composition cannot contain negative value at index " +
                QString::number(static_cast<int>(i))
                );
        }

        composition[static_cast<std::size_t>(component)] = value;
    }

    return composition;
}

std::array<double, ComponentCount> readComponentMassesKgArray(
    const QJsonArray& array,
    const std::vector<Component>& components
    )
{
    requireArraySize(
        array,
        components.size(),
        "initial_state.component_masses_kg"
        );

    std::array<double, ComponentCount> masses{};

    for (std::size_t i = 0; i < components.size(); ++i) {
        const Component component = components[i];

        const double value = readArrayNumber(
            array,
            static_cast<int>(i),
            "initial_state.component_masses_kg"
            );

        if (value < 0.0) {
            throw configError(
                "Component mass cannot be negative at index " +
                QString::number(static_cast<int>(i))
                );
        }

        masses[static_cast<std::size_t>(component)] = value;
    }

    return masses;
}

void validateCurrentComponentSet(const std::vector<Component>& components)
{
    if (components.size() != ComponentCount) {
        throw std::runtime_error(
            "At stage 1 config must contain exactly 3 components: C2H6, C5H12, H2O"
            );
    }

    for (std::size_t i = 0; i < ComponentCount; ++i) {
        const Component expected = static_cast<Component>(i);

        if (components[i] != expected) {
            throw std::runtime_error(
                "At stage 1 component order must be exactly: C2H6, C5H12, H2O"
                );
        }
    }
}

} // namespace

SimulationConfig SimulationConfigLoader::load(const std::string& path)
{
    QFile file(QString::fromStdString(path));

    if (!file.open(QIODevice::ReadOnly)) {
        const QString currentDir = QDir::currentPath();
        const QString absolutePath = QFileInfo(file).absoluteFilePath();

        throw std::runtime_error(
            "Cannot open simulation config file: " + path +
            "\nCurrent working directory: " + currentDir.toStdString() +
            "\nResolved config path: " + absolutePath.toStdString()
            );
    }

    const QByteArray data = file.readAll();

    QJsonParseError parseError;
    const QJsonDocument document = QJsonDocument::fromJson(data, &parseError);

    if (parseError.error != QJsonParseError::NoError) {
        throw std::runtime_error(
            "JSON parse error in " + path + ": " +
            parseError.errorString().toStdString()
            );
    }

    if (!document.isObject()) {
        throw std::runtime_error("Simulation config root must be a JSON object");
    }

    const QJsonObject root = document.object();

    SimulationConfig config;

    config.schemaVersion = readInt(root, "schema_version", 1);
    config.caseName = readString(root, "case_name", config.caseName);

    if (config.schemaVersion != 1) {
        throw std::runtime_error("Unsupported config schema_version");
    }

    const QJsonObject dataSources =
        requireObject(root, "data_sources");

    config.dataSources.componentPropertiesCsv =
        readString(
            dataSources,
            "component_properties_csv",
            "data/chemsep_components_properties.csv"
            );

    const QJsonArray componentsArray = requireArray(root, "components");

    config.componentKeys.clear();
    config.components.clear();

    for (const QJsonValue& value : componentsArray) {
        if (!value.isString()) {
            throw std::runtime_error("Every component name must be a string");
        }

        const QString componentKey = value.toString();

        config.componentKeys.push_back(
            componentKey.toStdString()
            );

        config.components.push_back(
            componentFromString(componentKey)
            );
    }

    validateCurrentComponentSet(config.components);

    const QJsonObject thermo = requireObject(root, "thermo");
    config.thermo.eos = readString(thermo, "eos", "PengRobinson");
    config.thermo.flashMode =
        flashModeFromString(QString::fromStdString(
            readString(thermo, "flash_mode", "RachfordRice")
            ));
    config.thermo.flashTolerance =
        readDouble(thermo, "flash_tolerance", 1.0e-8);
    config.thermo.flashMaxIterations =
        readInt(thermo, "flash_max_iterations", 100);

    if (config.thermo.eos != "PengRobinson") {
        throw std::runtime_error("At stage 1 only PengRobinson EOS is supported");
    }

    const QJsonObject feed = requireObject(root, "feed");
    config.feed.pressureBar = readDouble(feed, "pressure_bar", 2.5);
    config.feed.temperatureC = readDouble(feed, "temperature_C", 25.0);
    config.feed.flowKmolPerS = readDouble(feed, "flow_kmol_per_s", 0.001);
    config.feed.compositionBasis =
        readString(feed, "composition_basis", "mole_fraction");

    if (config.feed.compositionBasis != "mole_fraction") {
        throw std::runtime_error("At stage 1 only mole_fraction feed composition is supported");
    }

    config.feed.composition =
        readCompositionArray(
            requireArray(feed, "composition"),
            config.components
            );

    const QJsonObject reactor = requireObject(root, "reactor");
    config.reactor.type = readString(reactor, "type", "well_stirred");
    config.reactor.volumeM3 = readDouble(reactor, "volume_m3", 1.0);

    if (config.reactor.type != "well_stirred") {
        throw std::runtime_error("At stage 1 only well_stirred reactor type is supported");
    }

    const QJsonObject outlet = requireObject(reactor, "outlet");
    config.reactor.outlet.valveCoefficientKmolPerSBar =
        readDouble(outlet, "valve_coefficient_kmol_per_s_bar", 0.1);
    config.reactor.outlet.referencePressureBar =
        readDouble(outlet, "reference_pressure_bar", 2.5);
    config.reactor.outlet.allowReverseFlow =
        readBool(outlet, "allow_reverse_flow", false);

    const QJsonObject heatTransfer = requireObject(reactor, "heat_transfer");
    config.reactor.heatTransfer.areaM2 =
        readDouble(heatTransfer, "area_m2", 1.0);
    config.reactor.heatTransfer.coefficientWPerM2K =
        readDouble(heatTransfer, "coefficient_W_per_m2_K", 500.0);
    config.reactor.heatTransfer.jacketTemperatureC =
        readDouble(heatTransfer, "jacket_temperature_C", 110.0);

    const QJsonObject initialState = requireObject(root, "initial_state");
    config.initialState.temperatureC =
        readDouble(initialState, "temperature_C", 100.0);
    config.initialState.pressureBar =
        readDouble(initialState, "pressure_bar", 2.5);
    config.initialState.componentMassesKg =
        readComponentMassesKgArray(
            requireArray(initialState, "component_masses_kg"),
            config.components
            );

    const QJsonObject pressureInitialization =
        requireObject(initialState, "pressure_initialization");
    config.initialState.pressureInitialization.pressureMinPa =
        readDouble(pressureInitialization, "pressure_min_Pa", 1.0e4);
    config.initialState.pressureInitialization.pressureMaxPa =
        readDouble(pressureInitialization, "pressure_max_Pa", 1.0e7);
    config.initialState.pressureInitialization.volumeToleranceM3 =
        readDouble(pressureInitialization, "volume_tolerance_m3", 1.0e-8);
    config.initialState.pressureInitialization.newtonMaxIterations =
        readInt(pressureInitialization, "newton_max_iterations", 20);
    config.initialState.pressureInitialization.fallbackToBisection =
        readBool(pressureInitialization, "fallback_to_bisection", true);

    const QJsonObject reactions = requireObject(root, "reactions");
    config.reactions.enabled = readBool(reactions, "enabled", true);
    config.reactions.forwardRateConstant1PerS =
        readDouble(reactions, "forward_rate_constant_1_per_s", 100.0);
    config.reactions.reverseRateConstant1PerS =
        readDouble(reactions, "reverse_rate_constant_1_per_s", 100.0);

    const QJsonObject solver = requireObject(root, "solver");
    config.solver.model = readString(solver, "model", "energy_dae");
    config.solver.integrator = readString(solver, "integrator", "RadauIIA3");
    config.solver.timeStepS = readDouble(solver, "time_step_s", 1.0);
    config.solver.stepCount = readInt(solver, "step_count", 30);
    config.solver.printEvery = readInt(solver, "print_every", 1);

    if (config.solver.model != "energy_dae") {
        throw std::runtime_error("At stage 1 only energy_dae solver model is supported");
    }

    if (config.solver.integrator != "RadauIIA3") {
        throw std::runtime_error("At stage 1 only RadauIIA3 integrator is supported");
    }

    const QJsonObject limits = requireObject(root, "limits");

    config.limits.minComponentMassKg =
        readDouble(limits, "min_component_mass_kg", 0.0);
    config.limits.maxComponentMassKg =
        readDouble(limits, "max_component_mass_kg", 1.0e6);

    config.limits.minTemperatureC =
        readDouble(limits, "min_temperature_C", -100.0);
    config.limits.maxTemperatureC =
        readDouble(limits, "max_temperature_C", 1000.0);

    config.limits.minPressureBar =
        readDouble(limits, "min_pressure_bar", 0.001);
    config.limits.maxPressureBar =
        readDouble(limits, "max_pressure_bar", 1000.0);

    config.limits.minFlowKmolPerS =
        readDouble(limits, "min_flow_kmol_per_s", 0.0);
    config.limits.maxFlowKmolPerS =
        readDouble(limits, "max_flow_kmol_per_s", 1000.0);

    config.limits.minReactorVolumeM3 =
        readDouble(limits, "min_reactor_volume_m3", 1.0e-9);
    config.limits.maxReactorVolumeM3 =
        readDouble(limits, "max_reactor_volume_m3", 1.0e6);

    config.limits.minTimeStepS =
        readDouble(limits, "min_time_step_s", 1.0e-12);
    config.limits.maxTimeStepS =
        readDouble(limits, "max_time_step_s", 1.0e6);

    config.limits.minStepCount =
        readInt(limits, "min_step_count", 1);
    config.limits.maxStepCount =
        readInt(limits, "max_step_count", 10000000);

    config.limits.compositionSumTolerance =
        readDouble(limits, "composition_sum_tolerance", 1.0e-8);

    config.limits.failOnNegativeInventory =
        readBool(limits, "fail_on_negative_inventory", true);
    config.limits.failOnInvalidComposition =
        readBool(limits, "fail_on_invalid_composition", true);

    const QJsonObject output = requireObject(root, "output");

    config.output.enabled =
        readBool(output, "enabled", true);

    config.output.directory =
        readString(output, "directory", "results");

    config.output.writeResultTable =
        readBool(output, "write_result_table", true);

    config.output.resultTableFile =
        readString(output, "result_table_file", "result_table.tsv");

    config.output.writeInitialState =
        readBool(output, "write_initial_state", true);

    config.output.initialStateFile =
        readString(output, "initial_state_file", "initial_state.txt");

    config.output.writePressureInitialization =
        readBool(output, "write_pressure_initialization", true);

    config.output.pressureInitializationFile =
        readString(
            output,
            "pressure_initialization_file",
            "pressure_initialization.txt"
            );

    config.output.writeFinalState =
        readBool(output, "write_final_state", true);

    config.output.finalStateFile =
        readString(output, "final_state_file", "final_state.txt");

    const QJsonObject diagnostics = requireObject(root, "diagnostics");
    config.diagnostics.printInitialState =
        readBool(diagnostics, "print_initial_state", true);
    config.diagnostics.printMaterials =
        readBool(diagnostics, "print_materials", false);
    config.diagnostics.runFlashDiagnostics =
        readBool(diagnostics, "run_flash_diagnostics", false);
    config.diagnostics.runReactorModelDiagnostic =
        readBool(diagnostics, "run_reactor_model_diagnostic", false);
    config.diagnostics.runPressureInitializationDiagnostic =
        readBool(diagnostics, "run_pressure_initialization_diagnostic", true);
    config.diagnostics.runMassBalanceDiagnostic =
        readBool(diagnostics, "run_mass_balance_diagnostic", false);
    config.diagnostics.runDaeResidualDiagnostic =
        readBool(diagnostics, "run_dae_residual_diagnostic", false);
    config.diagnostics.runEnergyBalanceDiagnostic =
        readBool(diagnostics, "run_energy_balance_diagnostic", false);
    config.diagnostics.runEnergyDaeResidualDiagnostic =
        readBool(diagnostics, "run_energy_dae_residual_diagnostic", false);
    config.diagnostics.runEnergyRadauSingleStepDiagnostic =
        readBool(diagnostics, "run_energy_radau_single_step_diagnostic", false);
    config.diagnostics.runEnergyRadauSimulationDiagnostic =
        readBool(diagnostics, "run_energy_radau_simulation_diagnostic", true);
    config.diagnostics.printFinalState =
        readBool(diagnostics, "print_final_state", true);

    SimulationConfigValidator::validate(config);

    return config;
}
