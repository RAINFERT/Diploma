#include "BinaryInteractionDatabase.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

std::string trimText(const std::string& text) {
    std::size_t first = 0;
    while (first < text.size() && std::isspace(static_cast<unsigned char>(text[first]))) {
        ++first;
    }

    std::size_t last = text.size();
    while (last > first && std::isspace(static_cast<unsigned char>(text[last - 1]))) {
        --last;
    }

    return text.substr(first, last - first);
}

} // namespace

BinaryInteractionDatabase::BinaryInteractionDatabase(
    const std::string& csvPath
    ) {
    std::ifstream file(csvPath);

    if (!file.is_open()) {
        throw std::runtime_error(
            "Cannot open binary interaction file: " + csvPath
            );
    }

    const char delimiter = ';';

    std::string headerLine;
    if (!std::getline(file, headerLine)) {
        throw std::runtime_error(
            "Binary interaction file is empty: " + csvPath
            );
    }

    std::vector<std::string> header = splitCsvLine(headerLine, delimiter);

    if (header.size() < 2) {
        throw std::runtime_error(
            "Binary interaction file has invalid header: " + csvPath
            );
    }

    // header[0] пустой, потому что первая ячейка в CSV — это угол таблицы.
    // Начиная с header[1] идут имена компонентов.
    for (std::size_t j = 1; j < header.size(); ++j) {
        const std::string name = trimText(header[j]);
        names_.push_back(name);
        indexByName_[normalizeName(name)] = names_.size() - 1;
    }

    std::string line;
    std::size_t rowNumber = 0;

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        ++rowNumber;

        std::vector<std::string> fields = splitCsvLine(line, delimiter);

        if (fields.size() != header.size()) {
            throw std::runtime_error(
                "Invalid number of columns in binary interaction row: " +
                std::to_string(rowNumber)
                );
        }

        const std::string rowName = trimText(fields[0]);

        std::vector<double> row;
        row.reserve(names_.size());

        for (std::size_t j = 1; j < fields.size(); ++j) {
            row.push_back(parseDouble(fields[j], rowName, names_[j - 1]));
        }

        kij_.push_back(row);
    }

    if (kij_.size() != names_.size()) {
        throw std::runtime_error(
            "Binary interaction matrix is not square"
            );
    }
}

BinaryInteractionMatrix BinaryInteractionDatabase::buildMatrix(
    const ComponentSet& componentSet
    ) const {
    const std::size_t n = componentSet.size();

    BinaryInteractionMatrix result(
        n,
        std::vector<double>(n, 0.0)
        );

    for (std::size_t i = 0; i < n; ++i) {
        const std::string rowKey = normalizeName(componentSet[i].name);

        const auto rowIt = indexByName_.find(rowKey);
        if (rowIt == indexByName_.end()) {
            throw std::runtime_error(
                "Component not found in binary interaction rows: " +
                componentSet[i].name
                );
        }

        for (std::size_t j = 0; j < n; ++j) {
            const std::string columnKey = normalizeName(componentSet[j].name);

            const auto colIt = indexByName_.find(columnKey);
            if (colIt == indexByName_.end()) {
                throw std::runtime_error(
                    "Component not found in binary interaction columns: " +
                    componentSet[j].name
                    );
            }

            result[i][j] = kij_[rowIt->second][colIt->second];
        }
    }

    return result;
}

std::string BinaryInteractionDatabase::normalizeName(
    const std::string& name
    ) {
    std::string result = trimText(name);

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

std::vector<std::string> BinaryInteractionDatabase::splitCsvLine(
    const std::string& line,
    char delimiter
    ) {
    std::vector<std::string> fields;
    std::string field;
    std::stringstream stream(line);

    while (std::getline(stream, field, delimiter)) {
        fields.push_back(trimText(field));
    }

    return fields;
}

double BinaryInteractionDatabase::parseDouble(
    const std::string& value,
    const std::string& rowName,
    const std::string& columnName
    ) {
    const std::string text = trimText(value);

    if (text.empty()) {
        throw std::runtime_error(
            "Empty kij value for pair: " + rowName + " / " + columnName
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
            "Cannot parse kij value '" + text +
            "' for pair: " + rowName + " / " + columnName
            );
    }
}
