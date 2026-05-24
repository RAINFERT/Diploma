#include "ComponentSet.h"

#include <stdexcept>

ComponentSet::ComponentSet(
    std::vector<std::string> componentKeys,
    std::vector<Component> legacyComponents,
    std::vector<ChemsepComponent> chemsepComponents,
    MaterialList materials,
    ComponentEnthalpyDataList enthalpyData
    )
    : componentKeys_(std::move(componentKeys)),
    legacyComponents_(std::move(legacyComponents)),
    chemsepComponents_(std::move(chemsepComponents)),
    materials_(std::move(materials)),
    enthalpyData_(std::move(enthalpyData))
{
    validate();
}

std::size_t ComponentSet::size() const
{
    return componentKeys_.size();
}

const std::string& ComponentSet::key(
    const std::size_t index
    ) const
{
    return componentKeys_.at(index);
}

Component ComponentSet::legacyComponent(
    const std::size_t index
    ) const
{
    return legacyComponents_.at(index);
}

const ChemsepComponent& ComponentSet::chemsepComponent(
    const std::size_t index
    ) const
{
    return chemsepComponents_.at(index);
}

const Material& ComponentSet::material(
    const std::size_t index
    ) const
{
    return materials_.at(index);
}

const ComponentEnthalpyData& ComponentSet::enthalpy(
    const std::size_t index
    ) const
{
    return enthalpyData_.at(index);
}

const std::vector<std::string>& ComponentSet::keys() const
{
    return componentKeys_;
}

const std::vector<Component>& ComponentSet::legacyComponents() const
{
    return legacyComponents_;
}

const std::vector<ChemsepComponent>& ComponentSet::chemsepComponents() const
{
    return chemsepComponents_;
}

const MaterialList& ComponentSet::materials() const
{
    return materials_;
}

const ComponentEnthalpyDataList& ComponentSet::enthalpyData() const
{
    return enthalpyData_;
}

void ComponentSet::validate() const
{
    if (componentKeys_.empty()) {
        throw std::runtime_error(
            "ComponentSet must contain at least one component"
            );
    }

    if (componentKeys_.size() != legacyComponents_.size()) {
        throw std::runtime_error(
            "ComponentSet validation failed: componentKeys size differs from legacyComponents size"
            );
    }

    if (componentKeys_.size() != chemsepComponents_.size()) {
        throw std::runtime_error(
            "ComponentSet validation failed: componentKeys size differs from chemsepComponents size"
            );
    }

    if (materials_.size() != componentKeys_.size()) {
        throw std::runtime_error(
            "ComponentSet validation failed: materials size differs from componentKeys size"
            );
    }

    if (enthalpyData_.size() != componentKeys_.size()) {
        throw std::runtime_error(
            "ComponentSet validation failed: enthalpyData size differs from componentKeys size"
            );
    }


}
