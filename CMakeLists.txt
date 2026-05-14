cmake_minimum_required(VERSION 3.16)

project(WellStirredReactor VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)

set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
        mainwindow.ui
)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(WellStirredReactor
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
        core/ReactorState.h core/ReactorState.cpp
        core/Component.h
        core/thermo/Material.h
        core/thermo/Material.cpp
        core/thermo/PengRobinsonEOS.h core/thermo/PengRobinsonEOS.cpp
        core/thermo/FlashCalculation.h core/thermo/FlashCalculation.cpp
        core/thermo/ThermoPackage.h core/thermo/ThermoPackage.cpp
        core/models/WellStirredReactorModel.h core/models/WellStirredReactorModel.cpp
        core/reactions/ReactionModel.h core/reactions/ReactionModel.cpp
        core/dae/ReducedReactorDae.h core/dae/ReducedReactorDae.cpp
        core/numerics/NewtonSolver.h core/numerics/NewtonSolver.cpp
        core/simulation/ReactorSimulation.h core/simulation/ReactorSimulation.cpp
        core/thermo/EnthalpyModel.h core/thermo/EnthalpyModel.cpp
        core/dae/EnergyReactorDae.h core/dae/EnergyReactorDae.cpp
        core/simulation/EnergyReactorSimulation.h core/simulation/EnergyReactorSimulation.cpp
        README
        .gitignore
    )
# Define target properties for Android with Qt 6 as:
#    set_property(TARGET WellStirredReactor APPEND PROPERTY QT_ANDROID_PACKAGE_SOURCE_DIR
#                 ${CMAKE_CURRENT_SOURCE_DIR}/android)
# For more information, see https://doc.qt.io/qt-6/qt-add-executable.html#target-creation
else()
    if(ANDROID)
        add_library(WellStirredReactor SHARED
            ${PROJECT_SOURCES}
        )
# Define properties for Android with Qt 5 after find_package() calls as:
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
    else()
        add_executable(WellStirredReactor
            ${PROJECT_SOURCES}
        )
    endif()
endif()

target_link_libraries(WellStirredReactor PRIVATE Qt${QT_VERSION_MAJOR}::Widgets)

# Qt for iOS sets MACOSX_BUNDLE_GUI_IDENTIFIER automatically since Qt 6.1.
# If you are developing for iOS or macOS you should consider setting an
# explicit, fixed bundle identifier manually though.
if(${QT_VERSION} VERSION_LESS 6.1.0)
  set(BUNDLE_ID_OPTION MACOSX_BUNDLE_GUI_IDENTIFIER com.example.WellStirredReactor)
endif()
set_target_properties(WellStirredReactor PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

include(GNUInstallDirs)
install(TARGETS WellStirredReactor
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(WellStirredReactor)
endif()
