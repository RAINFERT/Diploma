#pragma once

#include "Component.h"

#include <string>
#include <unordered_map>
#include <vector>

using BinaryInteractionMatrix = std::vector<std::vector<double>>;

class BinaryInteractionDatabase {
public:
    explicit BinaryInteractionDatabase(const std::string& csvPath);

    BinaryInteractionMatrix buildMatrix(const ComponentSet& componentSet) const;

private:
    std::vector<std::string> names_;
    std::unordered_map<std::string, std::size_t> indexByName_;
    BinaryInteractionMatrix kij_;

    static std::string normalizeName(const std::string& name);

    static std::vector<std::string> splitCsvLine(
        const std::string& line,
        char delimiter
        );

    static double parseDouble(
        const std::string& value,
        const std::string& rowName,
        const std::string& columnName
        );
};
