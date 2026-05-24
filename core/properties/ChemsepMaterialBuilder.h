#pragma once

#include "../Component.h"
#include "../thermo/Material.h"
#include "ChemsepComponent.h"

#include <string>
#include <vector>

Material makeMaterialFromChemsep(
    const std::string& componentKey,
    const ChemsepComponent& chemsepComponent
    );

MaterialList makeMaterialListFromChemsep(
    const std::vector<std::string>& componentKeys,
    const std::vector<const ChemsepComponent*>& chemsepComponents
    );
