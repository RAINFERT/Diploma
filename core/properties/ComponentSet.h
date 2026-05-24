#pragma once

#include "../Component.h"
#include "../thermo/Material.h"
#include "../thermo/EnthalpyModel.h"
#include "ChemsepComponent.h"

#include <cstddef>
#include <string>
#include <vector>

class ComponentSet {
public:
    ComponentSet(
        std::vector<std::string> componentKeys,
        std::vector<Component> legacyComponents,
        std::vector<ChemsepComponent> chemsepComponents,
        MaterialList materials,
        ComponentEnthalpyDataList enthalpyData
        );

    std::size_t size() const;

    const std::string& key(
        std::size_t index
        ) const;

    Component legacyComponent(
        std::size_t index
        ) const;

    const ChemsepComponent& chemsepComponent(
        std::size_t index
        ) const;

    const Material& material(
        std::size_t index
        ) const;

    const ComponentEnthalpyData& enthalpy(
        std::size_t index
        ) const;

    const std::vector<std::string>& keys() const;
    const std::vector<Component>& legacyComponents() const;
    const std::vector<ChemsepComponent>& chemsepComponents() const;

    const MaterialList& materials() const;
    const ComponentEnthalpyDataList& enthalpyData() const;

private:
    std::vector<std::string> componentKeys_;
    std::vector<Component> legacyComponents_;
    std::vector<ChemsepComponent> chemsepComponents_;

    MaterialList materials_;
    ComponentEnthalpyDataList enthalpyData_;

    void validate() const;
};
