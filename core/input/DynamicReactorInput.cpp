#include "DynamicReactorInput.h"

#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QString>

#include <stdexcept>

namespace {

QJsonObject requireObject(
    const QJsonObject& parent,
    const QString& key
    ) {
    if (!parent.contains(key) || !parent.value(key).isObject()) {
        throw std::runtime_error(
            "reactor_input.json: required object not found: " +
            key.toStdString()
            );
    }

    return parent.value(key).toObject();
}

double requireDouble(
    const QJsonObject& object,
    const QString& key
    ) {
    if (!object.contains(key) || !object.value(key).isDouble()) {
        throw std::runtime_error(
            "reactor_input.json: required numeric field not found: " +
            key.toStdString()
            );
    }

    return object.value(key).toDouble();
}

int requireInt(
    const QJsonObject& object,
    const QString& key
    ) {
    if (!object.contains(key) || !object.value(key).isDouble()) {
        throw std::runtime_error(
            "reactor_input.json: required integer field not found: " +
            key.toStdString()
            );
    }

    return object.value(key).toInt();
}

bool requireBool(
    const QJsonObject& object,
    const QString& key
    ) {
    if (!object.contains(key) || !object.value(key).isBool()) {
        throw std::runtime_error(
            "reactor_input.json: required boolean field not found: " +
            key.toStdString()
            );
    }

    return object.value(key).toBool();
}

std::string requireString(
    const QJsonObject& object,
    const QString& key
    ) {
    if (!object.contains(key) || !object.value(key).isString()) {
        throw std::runtime_error(
            "reactor_input.json: required string field not found: " +
            key.toStdString()
            );
    }

    return object.value(key).toString().toStdString();
}

std::vector<std::string> requireStringArray(
    const QJsonObject& object,
    const QString& key
    ) {
    if (!object.contains(key) || !object.value(key).isArray()) {
        throw std::runtime_error(
            "reactor_input.json: required string array not found: " +
            key.toStdString()
            );
    }

    const QJsonArray array = object.value(key).toArray();

    std::vector<std::string> result;
    result.reserve(static_cast<std::size_t>(array.size()));

    for (const QJsonValue& value : array) {
        if (!value.isString()) {
            throw std::runtime_error(
                "reactor_input.json: non-string value in array: " +
                key.toStdString()
                );
        }

        result.push_back(
            value.toString().toStdString()
            );
    }

    return result;
}

DynamicComposition requireDoubleArray(
    const QJsonObject& object,
    const QString& key
    ) {
    if (!object.contains(key) || !object.value(key).isArray()) {
        throw std::runtime_error(
            "reactor_input.json: required numeric array not found: " +
            key.toStdString()
            );
    }

    const QJsonArray array = object.value(key).toArray();

    DynamicComposition result;
    result.reserve(static_cast<std::size_t>(array.size()));

    for (const QJsonValue& value : array) {
        if (!value.isDouble()) {
            throw std::runtime_error(
                "reactor_input.json: non-numeric value in array: " +
                key.toStdString()
                );
        }

        result.push_back(
            value.toDouble()
            );
    }

    return result;
}

} // namespace

DynamicReactorInput DynamicReactorInputReader::readJson(
    const std::string& filePath
    ) {
    QFile file(QString::fromStdString(filePath));

    if (!file.open(QIODevice::ReadOnly)) {
        throw std::runtime_error(
            "Cannot open reactor input file: " + filePath
            );
    }

    const QByteArray data = file.readAll();

    QJsonParseError parseError;
    const QJsonDocument document =
        QJsonDocument::fromJson(
            data,
            &parseError
            );

    if (parseError.error != QJsonParseError::NoError) {
        throw std::runtime_error(
            "reactor_input.json parse error: " +
            parseError.errorString().toStdString()
            );
    }

    if (!document.isObject()) {
        throw std::runtime_error(
            "reactor_input.json root must be an object"
            );
    }

    const QJsonObject root =
        document.object();

    DynamicReactorInput input;

    input.componentDatabaseCsv =
        requireString(
            root,
            "component_database_csv"
            );

    input.binaryInteractionCsv =
        requireString(
            root,
            "binary_interaction_csv"
            );

    input.componentNames =
        requireStringArray(
            root,
            "components"
            );

    const QJsonObject initial =
        requireObject(
            root,
            "initial"
            );

    input.initial.pressurePa =
        requireDouble(
            initial,
            "pressure_Pa"
            );

    input.initial.temperatureK =
        requireDouble(
            initial,
            "temperature_K"
            );

    input.initial.totalMolesKmol =
        requireDouble(
            initial,
            "total_moles_kmol"
            );

    input.initial.composition =
        requireDoubleArray(
            initial,
            "composition"
            );

    const QJsonObject feed =
        requireObject(
            root,
            "feed"
            );

    input.reactorParameters.inletMolarFlowKmolPerS =
        requireDouble(
            feed,
            "molar_flow_kmol_per_s"
            );

    input.reactorParameters.inletPressurePa =
        requireDouble(
            feed,
            "pressure_Pa"
            );

    input.reactorParameters.inletTemperatureK =
        requireDouble(
            feed,
            "temperature_K"
            );

    input.reactorParameters.inletComposition =
        requireDoubleArray(
            feed,
            "composition"
            );

    const QJsonObject valve =
        requireObject(
            root,
            "valve"
            );

    input.reactorParameters.outletValveCoefficientKmolPerSPa =
        requireDouble(
            valve,
            "outlet_valve_coefficient_kmol_per_s_Pa"
            );

    input.reactorParameters.allowReverseOutletFlow =
        requireBool(
            valve,
            "allow_reverse_outlet_flow"
            );

    const QJsonObject heatTransfer =
        requireObject(
            root,
            "heat_transfer"
            );

    input.reactorParameters.heatTransferAreaM2 =
        requireDouble(
            heatTransfer,
            "area_m2"
            );

    input.reactorParameters.heatTransferCoefficientWPerM2K =
        requireDouble(
            heatTransfer,
            "coefficient_W_per_m2_K"
            );

    input.reactorParameters.jacketTemperatureK =
        requireDouble(
            heatTransfer,
            "jacket_temperature_K"
            );

    const QJsonObject flash =
        requireObject(
            root,
            "flash"
            );

    input.flashOptions.initialLiquidSplit =
        requireDouble(
            flash,
            "initial_liquid_split"
            );

    input.flashOptions.pseudoTimeStep =
        requireDouble(
            flash,
            "pseudo_time_step"
            );

    input.flashOptions.tolerance =
        requireDouble(
            flash,
            "tolerance"
            );

    input.flashOptions.maxIterations =
        requireInt(
            flash,
            "max_iterations"
            );

    const QJsonObject radau =
        requireObject(
            root,
            "radau"
            );

    input.radauOptions.residualTolerance =
        requireDouble(
            radau,
            "residual_tolerance"
            );

    input.radauOptions.stepTolerance =
        requireDouble(
            radau,
            "step_tolerance"
            );

    input.radauOptions.maxNewtonIterations =
        requireInt(
            radau,
            "max_newton_iterations"
            );

    const QJsonObject simulation =
        requireObject(
            root,
            "simulation"
            );

    input.simulationOptions.timeStepS =
        requireDouble(
            simulation,
            "time_step_s"
            );

    input.simulationOptions.stepCount =
        requireInt(
            simulation,
            "step_count"
            );

    input.simulationOptions.radauOptions =
        input.radauOptions;

    return input;
}
