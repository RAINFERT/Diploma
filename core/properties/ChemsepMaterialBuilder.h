#pragma once

#include "../Component.h"
#include "../thermo/Material.h"
#include "ChemsepComponent.h"

#include <string>
#include <vector>

Material makeMaterialFromChemsep(
    Component component,
    const ChemsepComponent& chemsepComponent
    );

MaterialList makeMaterialListFromChemsep(
    const std::vector<Component>& components,
    const std::vector<const ChemsepComponent*>& chemsepComponents
    );
