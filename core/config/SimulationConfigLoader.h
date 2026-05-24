#pragma once

#include "SimulationConfig.h"

#include <string>

class SimulationConfigLoader {
public:
    static SimulationConfig load(const std::string& path);
};
