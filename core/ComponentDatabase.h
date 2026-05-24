#pragma once

#include "Component.h"

#include <string>
#include <unordered_map>
#include <vector>

class ComponentDatabase {
public:
    explicit ComponentDatabase(const std::string& csvPath);

    ComponentSet select(const std::vector<std::string>& componentNames) const;

    bool contains(const std::string& componentName) const;

private:
    std::unordered_map<std::string, ComponentProperties> componentsByName_;

    static std::string normalizeName(const std::string& name);

    static std::vector<std::string> splitCsvLine(
        const std::string& line,
        char delimiter
        );

    static double parseRequiredDouble(
        const std::string& value,
        const std::string& columnName,
        const std::string& componentName
        );

    static std::size_t findColumn(
        const std::vector<std::string>& header,
        const std::string& columnName
        );

    static int parseOptionalInt(
        const std::string& value,
        int defaultValue
        );

    static double parseOptionalDouble(
        const std::string& value,
        double defaultValue
        );

    static bool tryFindColumn(
        const std::vector<std::string>& header,
        const std::string& columnName,
        std::size_t& index
        );

    static bool tryFindAnyColumn(
        const std::vector<std::string>& header,
        const std::vector<std::string>& columnNames,
        std::size_t& index
        );
};
