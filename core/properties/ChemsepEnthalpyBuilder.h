#pragma once

#include "../Component.h"
#include "../thermo/EnthalpyModel.h"
#include "ChemsepComponent.h"

#include <vector>

ComponentEnthalpyData makeEnthalpyDataFromChemsep(
    Component component,
    const ChemsepComponent& chemsepComponent
    );

ComponentEnthalpyDataList makeEnthalpyDataListFromChemsep(
    const std::vector<Component>& components,
    const std::vector<const ChemsepComponent*>& chemsepComponents
    );
