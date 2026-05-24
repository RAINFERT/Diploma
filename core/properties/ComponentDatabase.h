#pragma once

#include "ChemsepComponent.h"

#include <string>
#include <unordered_map>
#include <vector>

class ComponentDatabase {
public:
    void addComponent(
        const ChemsepComponent& component
        );

    const ChemsepComponent& resolve(
        const std::string& key
        ) const;

    std::vector<const ChemsepComponent*> resolveMany(
        const std::vector<std::string>& keys
        ) const;

    std::size_t size() const;

private:
    std::vector<ChemsepComponent> components_;
    std::unordered_map<std::string, std::size_t> indexByKey_;

    void indexComponent(
        const ChemsepComponent& component,
        std::size_t index
        );

    void addIndex(
        const std::string& key,
        std::size_t index
        );

    static std::string normalizeKey(
        const std::string& key
        );

    static std::string aliasToChemsepName(
        const std::string& key
        );
};
