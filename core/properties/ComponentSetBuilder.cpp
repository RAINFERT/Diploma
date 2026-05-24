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
    if (config.componentKeys.empty()) {
        throw std::runtime_error(
            "Cannot build ComponentSet: empty component key list"
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
            config.componentKeys,
            resolvedComponents
            );

    const ComponentEnthalpyDataList enthalpyData =
        makeEnthalpyDataListFromChemsep(
            config.componentKeys,
            resolvedComponents
            );

    return ComponentSet(
        config.componentKeys,
        chemsepComponents,
        materials,
        enthalpyData
        );
}
