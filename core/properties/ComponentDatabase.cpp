#include "ComponentDatabase.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

namespace {

std::string trim(
    const std::string& value
    )
{
    const auto begin = std::find_if(
        value.begin(),
        value.end(),
        [](unsigned char c) {
            return !std::isspace(c);
        }
        );

    const auto end = std::find_if(
                         value.rbegin(),
                         value.rend(),
                         [](unsigned char c) {
                             return !std::isspace(c);
                         }
                         ).base();

    if (begin >= end) {
        return {};
    }

    return std::string(begin, end);
}

} // namespace

void ComponentDatabase::addComponent(
    const ChemsepComponent& component
    )
{
    const std::size_t index = components_.size();

    components_.push_back(component);

    indexComponent(
        components_.back(),
        index
        );
}

const ChemsepComponent& ComponentDatabase::resolve(
    const std::string& key
    ) const
{
    const std::string directKey =
        normalizeKey(key);

    const auto directIt =
        indexByKey_.find(directKey);

    if (directIt != indexByKey_.end()) {
        return components_.at(directIt->second);
    }

    const std::string aliasKey =
        normalizeKey(aliasToChemsepName(key));

    const auto aliasIt =
        indexByKey_.find(aliasKey);

    if (aliasIt != indexByKey_.end()) {
        return components_.at(aliasIt->second);
    }

    std::ostringstream message;
    message
        << "Cannot resolve component in ChemSep database: "
        << key
        << ". Use exact component_name, formula_chemsep, cas_number, "
        << "chemsep_index, or supported alias.";

    throw std::runtime_error(message.str());
}

std::vector<const ChemsepComponent*> ComponentDatabase::resolveMany(
    const std::vector<std::string>& keys
    ) const
{
    std::vector<const ChemsepComponent*> result;
    result.reserve(keys.size());

    for (const std::string& key : keys) {
        result.push_back(
            &resolve(key)
            );
    }

    return result;
}

std::size_t ComponentDatabase::size() const
{
    return components_.size();
}

void ComponentDatabase::indexComponent(
    const ChemsepComponent& component,
    const std::size_t index
    )
{
    addIndex(component.field("component_name"), index);
    addIndex(component.field("formula_chemsep"), index);
    addIndex(component.field("cas_number"), index);
    addIndex(component.field("chemsep_index"), index);
    addIndex(component.field("chemsep_record_id"), index);
}

void ComponentDatabase::addIndex(
    const std::string& key,
    const std::size_t index
    )
{
    const std::string normalized =
        normalizeKey(key);

    if (normalized.empty()) {
        return;
    }

    if (indexByKey_.find(normalized) == indexByKey_.end()) {
        indexByKey_[normalized] = index;
    }
}

std::string ComponentDatabase::normalizeKey(
    const std::string& key
    )
{
    std::string normalized =
        trim(key);

    std::transform(
        normalized.begin(),
        normalized.end(),
        normalized.begin(),
        [](unsigned char c) {
            return static_cast<char>(
                std::tolower(c)
                );
        }
        );

    return normalized;
}

std::string ComponentDatabase::aliasToChemsepName(
    const std::string& key
    )
{
    const std::string normalized =
        normalizeKey(key);

    if (normalized == "c2h6") {
        return "Ethane";
    }

    if (normalized == "c5h12") {
        return "N-pentane";
    }

    if (normalized == "h2o") {
        return "Water";
    }

    return key;
}
