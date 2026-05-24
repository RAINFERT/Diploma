#include "ComponentDatabase.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {

std::string trim(const std::string& text) {
    std::size_t first = 0;

    // Убираем UTF-8 BOM, если он есть в начале файла.
    // BOM выглядит как три байта: EF BB BF.
    if (text.size() >= 3
        && static_cast<unsigned char>(text[0]) == 0xEF
        && static_cast<unsigned char>(text[1]) == 0xBB
        && static_cast<unsigned char>(text[2]) == 0xBF) {
        first = 3;
    }

    while (first < text.size()
           && std::isspace(static_cast<unsigned char>(text[first]))) {
        ++first;
    }

    std::size_t last = text.size();

    while (last > first
           && std::isspace(static_cast<unsigned char>(text[last - 1]))) {
        --last;
    }

    return text.substr(first, last - first);
}

} // namespace

ComponentDatabase::ComponentDatabase(const std::string& csvPath) {
    std::ifstream file(csvPath);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open component database file: " + csvPath);
    }

    std::string headerLine;
    if (!std::getline(file, headerLine)) {
        throw std::runtime_error("Component database file is empty: " + csvPath);
    }

    const char delimiter = ',';
    const std::vector<std::string> header = splitCsvLine(headerLine, delimiter);

    const std::size_t idxChemsep = findColumn(header, "chemsep_index");
    const std::size_t idxName    = findColumn(header, "component_name");
    const std::size_t idxFormula = findColumn(header, "formula_chemsep");
    const std::size_t idxCas     = findColumn(header, "cas_number");

    const std::size_t idxMW      = findColumn(header, "MW_kg_per_kmol");
    const std::size_t idxTc      = findColumn(header, "Tc_K");
    const std::size_t idxPc      = findColumn(header, "Pc_Pa");
    const std::size_t idxVc      = findColumn(header, "Vc_m3_per_kmol");
    const std::size_t idxZc      = findColumn(header, "Zc");
    const std::size_t idxOmega   = findColumn(header, "omega");
    const std::size_t idxTb      = findColumn(header, "Tb_K");

    std::size_t idxCpEqno = 0;
    std::size_t idxCpA = 0;
    std::size_t idxCpB = 0;
    std::size_t idxCpC = 0;
    std::size_t idxCpD = 0;
    std::size_t idxCpE = 0;
    std::size_t idxCpTmin = 0;
    std::size_t idxCpTmax = 0;

    const bool hasCpEqno =
        tryFindAnyColumn(
            header,
            {
                "ideal_gas_cp_rpp_eqno",
                "ideal_gas_cp_eqno"
            },
            idxCpEqno
            );

    const bool hasCpA =
        tryFindAnyColumn(
            header,
            {
                "ideal_gas_cp_rpp_A",
                "ideal_gas_cp_A"
            },
            idxCpA
            );

    const bool hasCpB =
        tryFindAnyColumn(
            header,
            {
                "ideal_gas_cp_rpp_B",
                "ideal_gas_cp_B"
            },
            idxCpB
            );

    const bool hasCpC =
        tryFindAnyColumn(
            header,
            {
                "ideal_gas_cp_rpp_C",
                "ideal_gas_cp_C"
            },
            idxCpC
            );

    const bool hasCpD =
        tryFindAnyColumn(
            header,
            {
                "ideal_gas_cp_rpp_D",
                "ideal_gas_cp_D"
            },
            idxCpD
            );

    const bool hasCpE =
        tryFindAnyColumn(
            header,
            {
                "ideal_gas_cp_rpp_E",
                "ideal_gas_cp_E"
            },
            idxCpE
            );

    const bool hasCpTmin =
        tryFindAnyColumn(
            header,
            {
                "ideal_gas_cp_rpp_Tmin_K",
                "ideal_gas_cp_Tmin_K"
            },
            idxCpTmin
            );

    const bool hasCpTmax =
        tryFindAnyColumn(
            header,
            {
                "ideal_gas_cp_rpp_Tmax_K",
                "ideal_gas_cp_Tmax_K"
            },
            idxCpTmax
            );

    std::size_t idxLiquidCpEqno = 0;
    std::size_t idxLiquidCpA = 0;
    std::size_t idxLiquidCpB = 0;
    std::size_t idxLiquidCpC = 0;
    std::size_t idxLiquidCpD = 0;
    std::size_t idxLiquidCpE = 0;
    std::size_t idxLiquidCpTmin = 0;
    std::size_t idxLiquidCpTmax = 0;

    const bool hasLiquidCpEqno =
        tryFindColumn(header, "liquid_cp_eqno", idxLiquidCpEqno);

    const bool hasLiquidCpA =
        tryFindColumn(header, "liquid_cp_A", idxLiquidCpA);

    const bool hasLiquidCpB =
        tryFindColumn(header, "liquid_cp_B", idxLiquidCpB);

    const bool hasLiquidCpC =
        tryFindColumn(header, "liquid_cp_C", idxLiquidCpC);

    const bool hasLiquidCpD =
        tryFindColumn(header, "liquid_cp_D", idxLiquidCpD);

    const bool hasLiquidCpE =
        tryFindColumn(header, "liquid_cp_E", idxLiquidCpE);

    const bool hasLiquidCpTmin =
        tryFindColumn(header, "liquid_cp_Tmin_K", idxLiquidCpTmin);

    const bool hasLiquidCpTmax =
        tryFindColumn(header, "liquid_cp_Tmax_K", idxLiquidCpTmax);

    if (!hasLiquidCpA || !hasLiquidCpB || !hasLiquidCpC || !hasLiquidCpD) {
        std::cerr
            << "Warning: liquid Cp columns were not fully found. "
            << "DynamicEnthalpyModel may use fallback liquid Cp."
            << std::endl;
    }

    if (!hasCpA || !hasCpB || !hasCpC || !hasCpD) {
        std::cerr
            << "Warning: ideal gas Cp polynomial columns were not fully found. "
            << "DynamicEnthalpyModel may use fallback Cp."
            << std::endl;
    }

    std::string line;
    std::size_t lineNumber = 1;

    while (std::getline(file, line)) {
        ++lineNumber;

        if (line.empty()) {
            continue;
        }

        const std::vector<std::string> fields = splitCsvLine(line, delimiter);

        if (fields.size() < header.size()) {
            std::cerr << "Warning: skipped incomplete CSV line "
                      << lineNumber << std::endl;
            continue;
        }

        try {
            ComponentProperties component;

            component.chemsepIndex = static_cast<std::size_t>(
                parseRequiredDouble(fields[idxChemsep], "chemsep_index", fields[idxName])
                );

            component.name = trim(fields[idxName]);
            component.formula = trim(fields[idxFormula]);
            component.casNumber = trim(fields[idxCas]);

            component.molarMassKgPerKmol = parseRequiredDouble(
                fields[idxMW],
                "MW_kg_per_kmol",
                component.name
                );

            component.criticalTemperatureK = parseRequiredDouble(
                fields[idxTc],
                "Tc_K",
                component.name
                );

            component.criticalPressurePa = parseRequiredDouble(
                fields[idxPc],
                "Pc_Pa",
                component.name
                );

            component.criticalVolumeM3PerKmol = parseRequiredDouble(
                fields[idxVc],
                "Vc_m3_per_kmol",
                component.name
                );

            component.criticalCompressibility = parseRequiredDouble(
                fields[idxZc],
                "Zc",
                component.name
                );

            component.acentricFactor = parseRequiredDouble(
                fields[idxOmega],
                "omega",
                component.name
                );

            component.normalBoilingTemperatureK = parseRequiredDouble(
                fields[idxTb],
                "Tb_K",
                component.name
                );

            if (hasCpEqno) {
                component.idealGasCpEqno =
                    parseOptionalInt(
                        fields[idxCpEqno],
                        0
                        );
            }

            if (hasCpA) {
                component.idealGasCpA =
                    parseOptionalDouble(
                        fields[idxCpA],
                        0.0
                        );
            }

            if (hasCpB) {
                component.idealGasCpB =
                    parseOptionalDouble(
                        fields[idxCpB],
                        0.0
                        );
            }

            if (hasCpC) {
                component.idealGasCpC =
                    parseOptionalDouble(
                        fields[idxCpC],
                        0.0
                        );
            }

            if (hasCpD) {
                component.idealGasCpD =
                    parseOptionalDouble(
                        fields[idxCpD],
                        0.0
                        );
            }

            if (hasCpE) {
                component.idealGasCpE =
                    parseOptionalDouble(
                        fields[idxCpE],
                        0.0
                        );
            }

            if (hasCpTmin) {
                component.idealGasCpTminK =
                    parseOptionalDouble(
                        fields[idxCpTmin],
                        0.0
                        );
            }

            if (hasCpTmax) {
                component.idealGasCpTmaxK =
                    parseOptionalDouble(
                        fields[idxCpTmax],
                        0.0
                        );
            }

            if (hasLiquidCpEqno) {
                component.liquidCpEqno =
                    parseOptionalInt(
                        fields[idxLiquidCpEqno],
                        0
                        );
            }

            if (hasLiquidCpA) {
                component.liquidCpA =
                    parseOptionalDouble(
                        fields[idxLiquidCpA],
                        0.0
                        );
            }

            if (hasLiquidCpB) {
                component.liquidCpB =
                    parseOptionalDouble(
                        fields[idxLiquidCpB],
                        0.0
                        );
            }

            if (hasLiquidCpC) {
                component.liquidCpC =
                    parseOptionalDouble(
                        fields[idxLiquidCpC],
                        0.0
                        );
            }

            if (hasLiquidCpD) {
                component.liquidCpD =
                    parseOptionalDouble(
                        fields[idxLiquidCpD],
                        0.0
                        );
            }

            if (hasLiquidCpE) {
                component.liquidCpE =
                    parseOptionalDouble(
                        fields[idxLiquidCpE],
                        0.0
                        );
            }

            if (hasLiquidCpTmin) {
                component.liquidCpTminK =
                    parseOptionalDouble(
                        fields[idxLiquidCpTmin],
                        0.0
                        );
            }

            if (hasLiquidCpTmax) {
                component.liquidCpTmaxK =
                    parseOptionalDouble(
                        fields[idxLiquidCpTmax],
                        0.0
                        );
            }

            const std::string key = normalizeName(component.name);

            if (componentsByName_.count(key) > 0) {
                throw std::runtime_error(
                    "Duplicate component name in database: " + component.name
                    );
            }

            componentsByName_[key] = component;

        } catch (const std::exception& error) {
            std::cerr
                << "Warning: skipped component on CSV line "
                << lineNumber
                << ": "
                << error.what()
                << std::endl;

            continue;
        }
    }

    if (componentsByName_.empty()) {
        throw std::runtime_error("No components were loaded from: " + csvPath);
    }
}

ComponentSet ComponentDatabase::select(
    const std::vector<std::string>& componentNames
    ) const {
    if (componentNames.empty()) {
        throw std::runtime_error("Component list is empty");
    }

    ComponentSet result;

    for (const std::string& requestedName : componentNames) {
        const std::string key = normalizeName(requestedName);

        const auto it = componentsByName_.find(key);
        if (it == componentsByName_.end()) {
            throw std::runtime_error(
                "Component not found in ChemSep database: " + requestedName
                );
        }

        const std::size_t index = result.components.size();

        result.components.push_back(it->second);
        result.indexByName[key] = index;
    }

    return result;
}

bool ComponentDatabase::contains(const std::string& componentName) const {
    return componentsByName_.count(normalizeName(componentName)) > 0;
}

std::string ComponentDatabase::normalizeName(const std::string& name) {
    std::string result = trim(name);

    std::transform(
        result.begin(),
        result.end(),
        result.begin(),
        [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        }
        );

    return result;
}

std::vector<std::string> ComponentDatabase::splitCsvLine(
    const std::string& line,
    char delimiter
    ) {
    std::vector<std::string> fields;
    std::string field;

    bool insideQuotes = false;

    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

        if (c == '"') {
            // Две кавычки подряд внутри quoted field означают одну кавычку:
            // "hello ""world""" -> hello "world"
            if (insideQuotes && i + 1 < line.size() && line[i + 1] == '"') {
                field.push_back('"');
                ++i;
            } else {
                insideQuotes = !insideQuotes;
            }
        } else if (c == delimiter && !insideQuotes) {
            fields.push_back(trim(field));
            field.clear();
        } else {
            field.push_back(c);
        }
    }

    if (insideQuotes) {
        throw std::runtime_error("Unclosed quoted field in CSV line: " + line);
    }

    fields.push_back(trim(field));

    return fields;
}

double ComponentDatabase::parseRequiredDouble(
    const std::string& value,
    const std::string& columnName,
    const std::string& componentName
    ) {
    const std::string text = trim(value);

    if (text.empty()) {
        throw std::runtime_error(
            "Empty numeric value in column '" + columnName +
            "' for component '" + componentName + "'"
            );
    }

    try {
        std::size_t parsedCharacters = 0;
        const double number = std::stod(text, &parsedCharacters);

        if (parsedCharacters != text.size()) {
            throw std::runtime_error("not fully parsed");
        }

        return number;
    } catch (...) {
        throw std::runtime_error(
            "Cannot parse numeric value '" + text +
            "' in column '" + columnName +
            "' for component '" + componentName + "'"
            );
    }
}

std::size_t ComponentDatabase::findColumn(
    const std::vector<std::string>& header,
    const std::string& columnName
    ) {
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (trim(header[i]) == columnName) {
            return i;
        }
    }

    throw std::runtime_error("Required CSV column not found: " + columnName);
}

int ComponentDatabase::parseOptionalInt(
    const std::string& value,
    int defaultValue
    ) {
    const std::string text = trim(value);

    if (text.empty()) {
        return defaultValue;
    }

    try {
        std::size_t parsedCharacters = 0;
        const int number = std::stoi(text, &parsedCharacters);

        if (parsedCharacters != text.size()) {
            return defaultValue;
        }

        return number;
    } catch (...) {
        return defaultValue;
    }
}

double ComponentDatabase::parseOptionalDouble(
    const std::string& value,
    double defaultValue
    ) {
    const std::string text = trim(value);

    if (text.empty()) {
        return defaultValue;
    }

    try {
        std::size_t parsedCharacters = 0;
        const double number = std::stod(text, &parsedCharacters);

        if (parsedCharacters != text.size()) {
            return defaultValue;
        }

        return number;
    } catch (...) {
        return defaultValue;
    }
}

bool ComponentDatabase::tryFindColumn(
    const std::vector<std::string>& header,
    const std::string& columnName,
    std::size_t& index
    ) {
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (trim(header[i]) == columnName) {
            index = i;
            return true;
        }
    }

    return false;
}

bool ComponentDatabase::tryFindAnyColumn(
    const std::vector<std::string>& header,
    const std::vector<std::string>& columnNames,
    std::size_t& index
    ) {
    for (const std::string& columnName : columnNames) {
        if (tryFindColumn(header, columnName, index)) {
            return true;
        }
    }

    return false;
}
