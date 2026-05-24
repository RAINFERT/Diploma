#pragma once

#include "../config/SimulationConfig.h"
#include "ComponentDatabase.h"
#include "ComponentSet.h"

class ComponentSetBuilder {
public:
    static ComponentSet build(
        const SimulationConfig& config,
        const ComponentDatabase& database
        );
};
