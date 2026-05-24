#include "ComponentSetBuilder.h"

#include "ChemsepEnthalpyBuilder.h"
#include "ChemsepMaterialBuilder.h"

#include <stdexcept>
#include <vector>

ComponentSet ComponentSetBuilder::build(
    const SimulationConfig& config,
    const ComponentDatabase& database
    )
{
    if (config.componentKeys.size() != config.components.size()) {
        throw std::runtime_error(
            "Cannot build ComponentSet: componentKeys size differs from components size"
            );
    }

    const std::vector<const ChemsepComponent*> resolvedComponents =
        database.resolveMany(
            config.componentKeys
            );

    std::vector<ChemsepComponent> chemsepComponents;
    chemsepComponents.reserve(
        resolvedComponents.size()
        );

    for (const ChemsepComponent* component : resolvedComponents) {
        if (component == nullptr) {
            throw std::runtime_error(
                "Cannot build ComponentSet: null ChemSep component pointer"
                );
        }

        chemsepComponents.push_back(
            *component
            );
    }

    const MaterialList materials =
        makeMaterialListFromChemsep(
            config.components,
            resolvedComponents
            );

    const ComponentEnthalpyDataList enthalpyData =
        makeEnthalpyDataListFromChemsep(
            config.components,
            resolvedComponents
            );

    return ComponentSet(
        config.componentKeys,
        config.components,
        chemsepComponents,
        materials,
        enthalpyData
        );
}
