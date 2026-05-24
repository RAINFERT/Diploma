#pragma once

#include "../Component.h"
#include "../thermo/EnthalpyModel.h"
#include "ChemsepComponent.h"

#include <vector>

ComponentEnthalpyData makeEnthalpyDataFromChemsep(
    const std::string& componentKey,
    const ChemsepComponent& chemsepComponent
    );

ComponentEnthalpyDataList makeEnthalpyDataListFromChemsep(
    const std::vector<std::string>& componentKeys,
    const std::vector<const ChemsepComponent*>& chemsepComponents
    );
