#pragma once

#include "../thermo/ThermoPackage.h"

#include <string>
#include <vector>

struct DataSourcesConfig {
    std::string componentPropertiesCsv =
        "data/chemsep_components_properties.csv";
};

struct ThermoConfig {
    std::string eos = "PengRobinson";
    ThermoFlashMode flashMode = ThermoFlashMode::RachfordRice;
    double flashTolerance = 1.0e-8;
    int flashMaxIterations = 100;
};

struct FeedConfig {
    double pressureBar = 2.5;
    double temperatureC = 25.0;
    double flowKmolPerS = 0.001;
    std::string compositionBasis = "mole_fraction";

    // Ordered by SimulationConfig::componentKeys / components.
    std::vector<double> composition{0.5, 0.5, 0.0};
};

struct ReactorOutletConfig {
    double valveCoefficientKmolPerSBar = 0.1;
    double referencePressureBar = 2.5;
    bool allowReverseFlow = false;
};

struct ReactorHeatTransferConfig {
    double areaM2 = 1.0;
    double coefficientWPerM2K = 500.0;
    double jacketTemperatureC = 110.0;
};

struct ReactorConfig {
    std::string type = "well_stirred";
    double volumeM3 = 1.0;
    ReactorOutletConfig outlet;
    ReactorHeatTransferConfig heatTransfer;
};

struct PressureInitializationConfig {
    double pressureMinPa = 1.0e4;
    double pressureMaxPa = 1.0e7;
    double volumeToleranceM3 = 1.0e-8;
    int newtonMaxIterations = 20;
    bool fallbackToBisection = true;
};

struct InitialStateConfig {
    double temperatureC = 100.0;
    double pressureBar = 2.5;

    // Ordered by SimulationConfig::componentKeys / components.
    std::vector<double> componentMassesKg{0.0, 1.0, 500.0};

    PressureInitializationConfig pressureInitialization;
};

struct ReactionsConfig {
    bool enabled = true;
    double forwardRateConstant1PerS = 100.0;
    double reverseRateConstant1PerS = 100.0;
};

struct SolverConfig {
    std::string model = "energy_dae";
    std::string integrator = "RadauIIA3";
    double timeStepS = 1.0;
    int stepCount = 30;
    int printEvery = 1;
};

struct OutputConfig {
    bool enabled = true;
    std::string directory = "results";

    bool writeResultTable = true;
    std::string resultTableFile = "result_table.csv";

    bool writeInitialState = true;
    std::string initialStateFile = "initial_state.txt";

    bool writePressureInitialization = true;
    std::string pressureInitializationFile = "pressure_initialization.txt";

    bool writeFinalState = true;
    std::string finalStateFile = "final_state.txt";
};

struct DiagnosticsConfig {
    bool printInitialState = true;
    bool printMaterials = false;
    bool runFlashDiagnostics = false;
    bool runReactorModelDiagnostic = false;
    bool runPressureInitializationDiagnostic = true;
    bool runMassBalanceDiagnostic = false;
    bool runDaeResidualDiagnostic = false;
    bool runEnergyBalanceDiagnostic = false;
    bool runEnergyDaeResidualDiagnostic = false;
    bool runEnergyRadauSingleStepDiagnostic = false;
    bool runEnergyRadauSimulationDiagnostic = true;
    bool printFinalState = true;
};

struct LimitsConfig {
    double minComponentMassKg = 0.0;
    double maxComponentMassKg = 1.0e6;

    double minTemperatureC = -100.0;
    double maxTemperatureC = 1000.0;

    double minPressureBar = 0.001;
    double maxPressureBar = 1000.0;

    double minFlowKmolPerS = 0.0;
    double maxFlowKmolPerS = 1000.0;

    double minReactorVolumeM3 = 1.0e-9;
    double maxReactorVolumeM3 = 1.0e6;

    double minTimeStepS = 1.0e-12;
    double maxTimeStepS = 1.0e6;

    int minStepCount = 1;
    int maxStepCount = 10000000;

    double compositionSumTolerance = 1.0e-8;

    bool failOnNegativeInventory = true;
    bool failOnInvalidComposition = true;
};

struct SimulationConfig {
    int schemaVersion = 1;
    std::string caseName = "current_slave_v1_baseline";

    DataSourcesConfig dataSources;

    // Ordered component identifiers from input file.
    // Resolved through ChemSep ComponentDatabase.
    std::vector<std::string> componentKeys;

    ThermoConfig thermo;
    FeedConfig feed;
    ReactorConfig reactor;
    InitialStateConfig initialState;
    ReactionsConfig reactions;
    SolverConfig solver;
    LimitsConfig limits;
    OutputConfig output;
    DiagnosticsConfig diagnostics;
};
