#include "ComponentDatabaseLoader.h"

#include <cstdlib>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace {

std::vector<std::string> parseCsvLine(
    const std::string& line
    )
{
    std::vector<std::string> fields;

    std::string current;
    bool insideQuotes = false;

    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

        if (c == '"') {
            if (
                insideQuotes &&
                i + 1 < line.size() &&
                line[i + 1] == '"'
                ) {
                current.push_back('"');
                ++i;
            }
            else {
                insideQuotes = !insideQuotes;
            }
        }
        else if (c == ',' && !insideQuotes) {
            fields.push_back(current);
            current.clear();
        }
        else {
            current.push_back(c);
        }
    }

    fields.push_back(current);

    return fields;
}

std::string removeUtf8Bom(
    const std::string& value
    )
{
    if (
        value.size() >= 3 &&
        static_cast<unsigned char>(value[0]) == 0xEF &&
        static_cast<unsigned char>(value[1]) == 0xBB &&
        static_cast<unsigned char>(value[2]) == 0xBF
        ) {
        return value.substr(3);
    }

    return value;
}

std::string trim(
    const std::string& value
    )
{
    const std::string withoutBom =
        removeUtf8Bom(value);

    const auto first =
        withoutBom.find_first_not_of(" \t\r\n");

    if (first == std::string::npos) {
        return {};
    }

    const auto last =
        withoutBom.find_last_not_of(" \t\r\n");

    return withoutBom.substr(
        first,
        last - first + 1
        );
}

bool tryParseDouble(
    const std::string& text,
    double& value
    )
{
    const std::string cleaned =
        trim(text);

    if (cleaned.empty()) {
        return false;
    }

    char* endPointer = nullptr;

    value = std::strtod(
        cleaned.c_str(),
        &endPointer
        );

    if (endPointer == cleaned.c_str()) {
        return false;
    }

    while (*endPointer != '\0') {
        if (!std::isspace(
                static_cast<unsigned char>(*endPointer)
                )) {
            return false;
        }

        ++endPointer;
    }

    return std::isfinite(value);
}

std::unordered_map<std::string, std::string> makeRowMap(
    const std::vector<std::string>& headers,
    const std::vector<std::string>& values
    )
{
    std::unordered_map<std::string, std::string> row;

    for (std::size_t i = 0; i < headers.size(); ++i) {
        const std::string value =
            i < values.size()
                ? values[i]
                : std::string();

        row[headers[i]] = value;
    }

    return row;
}

std::vector<std::string> correlationNames()
{
    return {
        "antoine",
        "vapour_pressure",
        "heat_vaporization",
        "ideal_gas_cp",
        "ideal_gas_cp_rpp",
        "liquid_cp",
        "solid_cp",
        "liquid_density",
        "solid_density",
        "liquid_viscosity",
        "liquid_viscosity_rps",
        "vapour_viscosity",
        "liquid_thermal_conductivity",
        "vapour_thermal_conductivity",
        "surface_tension",
        "second_virial",
        "relative_static_permittivity"
    };
}

double optionalDouble(
    const std::unordered_map<std::string, std::string>& row,
    const std::string& key
    )
{
    const auto it = row.find(key);

    if (it == row.end()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double value = std::numeric_limits<double>::quiet_NaN();

    if (!tryParseDouble(it->second, value)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return value;
}

int optionalInt(
    const std::unordered_map<std::string, std::string>& row,
    const std::string& key
    )
{
    const double value =
        optionalDouble(row, key);

    if (!std::isfinite(value)) {
        return -1;
    }

    return static_cast<int>(value);
}

void addCorrelations(
    ChemsepComponent& component,
    const std::unordered_map<std::string, std::string>& row
    )
{
    for (const std::string& name : correlationNames()) {
        ChemsepCorrelation correlation;

        correlation.eqno =
            optionalInt(row, name + "_eqno");

        const auto unitIt =
            row.find(name + "_unit");

        if (unitIt != row.end()) {
            correlation.unit = unitIt->second;
        }

        correlation.A = optionalDouble(row, name + "_A");
        correlation.B = optionalDouble(row, name + "_B");
        correlation.C = optionalDouble(row, name + "_C");
        correlation.D = optionalDouble(row, name + "_D");
        correlation.E = optionalDouble(row, name + "_E");

        correlation.TminK =
            optionalDouble(row, name + "_Tmin_K");

        correlation.TmaxK =
            optionalDouble(row, name + "_Tmax_K");

        correlation.hasData =
            correlation.eqno >= 0;

        if (correlation.hasData) {
            component.setCorrelation(
                name,
                correlation
                );
        }
    }
}

ChemsepComponent makeComponent(
    const std::vector<std::string>& headers,
    const std::vector<std::string>& values
    )
{
    const auto row =
        makeRowMap(headers, values);

    ChemsepComponent component;

    for (const std::string& header : headers) {
        const auto it =
            row.find(header);

        const std::string value =
            it == row.end()
                ? std::string()
                : it->second;

        component.setField(
            header,
            value
            );

        double numericValue = 0.0;

        if (tryParseDouble(value, numericValue)) {
            component.setScalar(
                header,
                numericValue
                );
        }
    }

    addCorrelations(
        component,
        row
        );

    return component;
}

void requireHeaders(
    const std::vector<std::string>& headers
    )
{
    const std::unordered_set<std::string> required = {
        "chemsep_index",
        "component_name",
        "formula_chemsep",
        "cas_number",
        "MW_kg_per_kmol",
        "Tc_K",
        "Pc_Pa",
        "Vc_m3_per_kmol",
        "Zc",
        "omega",
        "h_form_ig_J_per_kmol"
    };

    std::unordered_set<std::string> actual;

    for (const std::string& header : headers) {
        actual.insert(header);
    }

    for (const std::string& header : required) {
        if (actual.find(header) == actual.end()) {
            throw std::runtime_error(
                "Required ChemSep column is missing: " +
                header
                );
        }
    }
}

} // namespace

ComponentDatabase ComponentDatabaseLoader::loadFromCsv(
    const std::string& path
    )
{
    std::ifstream file(path);

    if (!file.is_open()) {
        throw std::runtime_error(
            "Cannot open ChemSep component properties CSV: " +
            path
            );
    }

    std::string headerLine;

    if (!std::getline(file, headerLine)) {
        throw std::runtime_error(
            "ChemSep component properties CSV is empty: " +
            path
            );
    }

    std::vector<std::string> headers =
        parseCsvLine(headerLine);

    for (std::string& header : headers) {
        header = trim(header);
    }

    requireHeaders(headers);

    ComponentDatabase database;

    std::string line;
    std::size_t lineNumber = 1;

    while (std::getline(file, line)) {
        ++lineNumber;

        if (trim(line).empty()) {
            continue;
        }

        const std::vector<std::string> values =
            parseCsvLine(line);

        if (values.size() > headers.size()) {
            std::ostringstream message;
            message
                << "Too many fields in ChemSep CSV at line "
                << lineNumber
                << ". Expected at most "
                << headers.size()
                << ", got "
                << values.size();

            throw std::runtime_error(message.str());
        }

        ChemsepComponent component =
            makeComponent(headers, values);

        if (component.componentName().empty()) {
            std::ostringstream message;
            message
                << "ChemSep component_name is empty at line "
                << lineNumber;

            throw std::runtime_error(message.str());
        }

        database.addComponent(component);
    }

    if (database.size() == 0) {
        throw std::runtime_error(
            "ChemSep component database is empty after loading: " +
            path
            );
    }

    return database;
}
