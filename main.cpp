#include "mainwindow.h"

#include "core/ReactorState.h"
#include "core/models/WellStirredReactorModel.h"
#include "core/reactions/ReactionModel.h"
#include "core/thermo/ThermoPackage.h"
#include "core/dae/EnergyReactorDae.h"
#include "core/simulation/EnergyReactorSimulation.h"
#include "core/config/SimulationConfigLoader.h"
#include "core/properties/ComponentDatabaseLoader.h"
#include "core/properties/ComponentSetBuilder.h"

#include <QApplication>
#include <QDir>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

namespace
{
    WellStirredReactorParameters makeReactorParameters(
        const SimulationConfig& config
    )
    {
        WellStirredReactorParameters parameters;

        parameters.volumeM3 =
            config.reactor.volumeM3;

        parameters.inletFlowKmolPerS =
            config.feed.flowKmolPerS;

        parameters.inletPressureBar =
            config.feed.pressureBar;

        parameters.inletTemperatureC =
            config.feed.temperatureC;

        parameters.inletComposition =
            config.feed.composition;

        parameters.outletValveCoefficientKmolPerSBar =
            config.reactor.outlet.valveCoefficientKmolPerSBar;

        parameters.outletReferencePressureBar =
            config.reactor.outlet.referencePressureBar;

        parameters.allowReverseOutletFlow =
            config.reactor.outlet.allowReverseFlow;

        parameters.heatTransferAreaM2 =
            config.reactor.heatTransfer.areaM2;

        parameters.heatTransferCoefficientWPerM2K =
            config.reactor.heatTransfer.coefficientWPerM2K;

        parameters.jacketTemperatureC =
            config.reactor.heatTransfer.jacketTemperatureC;

        return parameters;
    }

    ReactorState makeInitialState(
        const SimulationConfig& config,
        const ComponentSet& componentSet
        )
    {
        std::vector<std::string> names;
        std::vector<double> molarMasses;

        names.reserve(componentSet.size());
        molarMasses.reserve(componentSet.size());

        for (std::size_t i = 0; i < componentSet.size(); ++i) {
            names.push_back(
                componentSet.key(i)
                );

            molarMasses.push_back(
                componentSet.material(i).molarMassKgPerKmol
                );
        }

        ReactorState state(
            names,
            molarMasses
            );

        if (config.initialState.componentMassesKg.size() != componentSet.size()) {
            throw std::runtime_error(
                "Initial component mass vector size differs from ComponentSet size"
                );
        }

        for (std::size_t i = 0; i < componentSet.size(); ++i) {
            state.setMassKg(
                i,
                config.initialState.componentMassesKg[i]
                );
        }

        state.setTemperatureC(
            config.initialState.temperatureC
            );

        state.setPressureBar(
            config.initialState.pressureBar
            );

        state.setEnergyJ(0.0);

        return state;
    }

    EnergyReactorSimulationOptions makeEnergySimulationOptions(
        const SimulationConfig& config
    )
    {
        EnergyReactorSimulationOptions options;

        options.timeStepS =
            config.solver.timeStepS;

        options.stepCount =
            config.solver.stepCount;

        return options;
    }

    void printSectionTitle(const std::string& title)
    {
        std::cout << "\n========================================\n";
        std::cout << title << "\n";
        std::cout << "========================================\n";
    }

    std::string pressureInitializationSummaryToString(
        const PressureInitializationResult& result
        )
    {
        std::ostringstream stream;

        stream << std::fixed << std::setprecision(8);
        stream << "converged  = "
               << (result.converged ? "true" : "false") << "\n";
        stream << "method     = " << result.method << "\n";
        stream << "iterations = " << result.iterations << "\n";
        stream << "P          = " << result.pressureBar << " bar\n";
        stream << "residual   = " << result.residualM3 << " m3\n";
        stream << "message    = " << result.message << "\n";

        return stream.str();
    }

    void writeTextFile(
        const std::string& path,
        const std::string& text
        )
    {
        std::ofstream file(path);

        if (!file.is_open()) {
            throw std::runtime_error(
                "Cannot open output file for writing: " + path
                );
        }

        file << text;
    }

    std::string makeOutputPath(
        const OutputConfig& output,
        const std::string& fileName
        )
    {
        QDir directory(
            QString::fromStdString(output.directory)
            );

        return directory.filePath(
                            QString::fromStdString(fileName)
                            ).toStdString();
    }

    void ensureOutputDirectory(
        const OutputConfig& output
        )
    {
        if (!output.enabled) {
            return;
        }

        QDir directory;

        const bool created =
            directory.mkpath(
                QString::fromStdString(output.directory)
                );

        if (!created) {
            throw std::runtime_error(
                "Cannot create output directory: " +
                output.directory
                );
        }
    }

    void printPressureInitializationSummary(
        const PressureInitializationResult& result
        )
    {
        printSectionTitle("Pressure initialization summary");

        std::cout << pressureInitializationSummaryToString(result);
    }

    void printSelectedChemsepComponents(
        const ComponentSet& componentSet
        )
    {
        printSectionTitle("Selected ChemSep components");

        std::cout
            << "index\tinput_key\tchemsep_name\tformula\tCAS\tMW\tTc_K\tPc_Pa\tomega\n";

        for (std::size_t i = 0; i < componentSet.size(); ++i) {
            const ChemsepComponent& component =
                componentSet.chemsepComponent(i);

            std::cout
                << i << "\t"
                << componentSet.key(i) << "\t"
                << component.componentName() << "\t"
                << component.formulaChemsep() << "\t"
                << component.casNumber() << "\t"
                << component.molarMassKgPerKmol() << "\t"
                << component.criticalTemperatureK() << "\t"
                << component.criticalPressurePa() << "\t"
                << component.acentricFactor()
                << "\n";
        }

        std::cout << std::endl;
    }

    void printThermoPackageMaterials(
        const ComponentSet& componentSet
        )
    {
        printSectionTitle("ThermoPackage materials");

        std::cout
            << "index\tname\tMW\tTc\tPc\tVc\tZc\tomega\n";

        for (std::size_t i = 0; i < componentSet.size(); ++i) {
            const Material& material =
                componentSet.material(i);

            std::cout
                << i << "\t"
                << material.name << "\t"
                << material.molarMassKgPerKmol << "\t"
                << material.criticalTemperatureK << "\t"
                << material.criticalPressurePa << "\t"
                << material.criticalVolumeM3PerKmol << "\t"
                << material.criticalCompressibility << "\t"
                << material.acentricFactor
                << "\n";
        }

        std::cout << std::endl;
    }

    void printThermoPackageEnthalpyData(
        const ComponentSet& componentSet
        )
    {
        printSectionTitle("ThermoPackage enthalpy data");

        std::cout
            << "index\tname\tHref_J_per_kmol\tCpA\tCpB\tCpC\tCpD\tCpE\n";

        for (std::size_t i = 0; i < componentSet.size(); ++i) {
            const ComponentEnthalpyData& item =
                componentSet.enthalpy(i);

            std::cout
                << i << "\t"
                << item.name << "\t"
                << item.referenceMolarEnthalpyJPerKmol << "\t"
                << item.cpA << "\t"
                << item.cpB << "\t"
                << item.cpC << "\t"
                << item.cpD << "\t"
                << item.cpE
                << "\n";
        }

        std::cout << std::endl;
    }
}


int main(int argc, char* argv[])
{
    try
    {
        const std::string configPath =
            (argc > 1)
                ? std::string(argv[1])
                : std::string("input/case_default.json");

        const SimulationConfig config =
            SimulationConfigLoader::load(configPath);

        ensureOutputDirectory(config.output);

        const ComponentDatabase componentDatabase =
            ComponentDatabaseLoader::loadFromCsv(
                config.dataSources.componentPropertiesCsv
                );

        const ComponentSet componentSet =
            ComponentSetBuilder::build(
                config,
                componentDatabase
                );

        std::cout << "Loaded ChemSep component database: "
                  << componentDatabase.size()
                  << " components\n";

        if (config.diagnostics.printMaterials) {
            printSelectedChemsepComponents(
                componentSet
                );
        }

        std::cout << "Loaded simulation config: "
                  << config.caseName << "\n";

        std::cout << "Config flash mode = "
                  << thermoFlashModeToString(config.thermo.flashMode)
                  << "\n";

        const ThermoFlashMode flashMode =
            config.thermo.flashMode;

        ThermoPackage thermo(
            componentSet.materials(),
            componentSet.enthalpyData(),
            flashMode
            );

        if (config.diagnostics.printMaterials) {
            printThermoPackageMaterials(
                componentSet
                );

            printThermoPackageEnthalpyData(
                componentSet
                );
        }

        std::cout << "Thermo flash mode = "
                  << thermoFlashModeToString(thermo.flashMode())
                  << "\n";

        ReactorState initialState =
            makeInitialState(
                config,
                componentSet
                );

        ReactionModel reactionModel(
            thermo.materials()
        );

        reactionModel.setEnabled(
            config.reactions.enabled
            );

        if (config.reactions.enabled)
        {
            reactionModel.setForwardRateConstant(
                config.reactions.forwardRateConstant1PerS
            );

            reactionModel.setReverseRateConstant(
                config.reactions.reverseRateConstant1PerS
            );
        }
        else
        {
            reactionModel.setForwardRateConstant(0.0);
            reactionModel.setReverseRateConstant(0.0);
        }

        const WellStirredReactorParameters parameters =
            makeReactorParameters(config);

        const WellStirredReactorModel reactorModel(
            parameters,
            thermo,
            reactionModel
        );

        if (config.diagnostics.printInitialState)
        {
            printSectionTitle("Initial reactor state");
            std::cout << initialState.toString() << std::endl;
        }

        if (
            config.output.enabled &&
            config.output.writeInitialState
            )
        {
            writeTextFile(
                makeOutputPath(
                    config.output,
                    config.output.initialStateFile
                    ),
                initialState.toString()
                );
        }

        const PressureInitializationConfig& pressureConfig =
            config.initialState.pressureInitialization;

        const PressureInitializationResult pressureInitialization =
            reactorModel.initializePressureForVolumeNewton(
                initialState,
                pressureConfig.pressureMinPa,
                pressureConfig.pressureMaxPa,
                pressureConfig.volumeToleranceM3,
                pressureConfig.newtonMaxIterations,
                pressureConfig.fallbackToBisection
            );

        if (config.diagnostics.runPressureInitializationDiagnostic)
        {
            printPressureInitializationSummary(
                pressureInitialization
            );
        }

        if (!pressureInitialization.converged)
        {
            std::cerr << "Pressure initialization failed. Stop simulation.\n";
            return 1;
        }

        if (
            config.output.enabled &&
            config.output.writePressureInitialization
            )
        {
            writeTextFile(
                makeOutputPath(
                    config.output,
                    config.output.pressureInitializationFile
                    ),
                pressureInitializationSummaryToString(
                    pressureInitialization
                    )
                );
        }

        const ReactorState initializedState =
            reactorModel.initializeEnergyFromEnthalpy(
                pressureInitialization.initializedState
            );

        EnergyReactorDae energyDae(
            reactorModel,
            initializedState
        );

        EnergyReactorSimulationRunner energySimulationRunner(
            energyDae,
            reactorModel
        );

        const EnergyReactorSimulationOptions energySimulationOptions =
            makeEnergySimulationOptions(config);

        const EnergyReactorSimulationResult energySimulationResult =
            energySimulationRunner.runRadauIIA3(
                initializedState,
                energySimulationOptions
            );

        const std::string resultTable =
            energySimulationRunner.resultTableToString(
                energySimulationResult,
                config.solver.printEvery
                );

        if (config.diagnostics.runEnergyRadauSimulationDiagnostic)
        {
            printSectionTitle(
                "Energy Radau IIA 3-stage DAE simulation diagnostic"
                );

            std::cout << "dt = "
                      << energySimulationOptions.timeStepS
                      << " s\n";

            std::cout << "steps = "
                      << energySimulationOptions.stepCount
                      << "\n\n";

            std::cout << resultTable << std::endl;
        }

        if (
            config.output.enabled &&
            config.output.writeResultTable
            )
        {
            writeTextFile(
                makeOutputPath(
                    config.output,
                    config.output.resultTableFile
                    ),
                resultTable
                );
        }

        if (!energySimulationResult.converged)
        {
            std::cerr << "Energy simulation failed: "
                      << energySimulationResult.message
                      << std::endl;

            return 1;
        }

        ReactorState finalState =
            initializedState;

        if (!energySimulationResult.points.empty())
        {
            finalState =
                energySimulationResult.finalState();
        }

        if (config.diagnostics.printFinalState)
        {
            printSectionTitle("Final reactor state");
            std::cout << finalState.toString() << std::endl;
        }

        if (
            config.output.enabled &&
            config.output.writeFinalState
            )
        {
            writeTextFile(
                makeOutputPath(
                    config.output,
                    config.output.finalStateFile
                    ),
                finalState.toString()
                );
        }

        QApplication app(argc, argv);

        MainWindow window;
        window.show();

        return app.exec();
    }
    catch (const std::exception& exception)
    {
        std::cerr << "\nFatal error:\n"
                  << exception.what()
                  << "\n";

        return 1;
    }
}
