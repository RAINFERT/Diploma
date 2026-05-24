#pragma once

#include "../Component.h"
#include "../dae/DynamicEnergyReactorRadau.h"
#include "../models/DynamicWellStirredReactorModel.h"
#include "../simulation/DynamicEnergyReactorSimulation.h"
#include "../thermo/DynamicGibbsFlashCalculation.h"

#include <string>
#include <vector>

struct DynamicInitialInput {
    double pressurePa = 101325.0;
    double temperatureK = 298.15;
    double totalMolesKmol = 1.0;
    DynamicComposition composition;
};

struct DynamicReactorInput {
    std::string componentDatabaseCsv;
    std::string binaryInteractionCsv;

    std::vector<std::string> componentNames;

    DynamicInitialInput initial;

    DynamicWellStirredReactorParameters reactorParameters;

    DynamicGibbsFlashOptions flashOptions;
    DynamicRadauOptions radauOptions;
    DynamicEnergyReactorSimulationOptions simulationOptions;
};

class DynamicReactorInputReader {
public:
    static DynamicReactorInput readJson(
        const std::string& filePath
        );
};
