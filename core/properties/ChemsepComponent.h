#pragma once

#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

struct ChemsepCorrelation {
    int eqno = -1;
    std::string unit;

    double A = std::numeric_limits<double>::quiet_NaN();
    double B = std::numeric_limits<double>::quiet_NaN();
    double C = std::numeric_limits<double>::quiet_NaN();
    double D = std::numeric_limits<double>::quiet_NaN();
    double E = std::numeric_limits<double>::quiet_NaN();

    double TminK = std::numeric_limits<double>::quiet_NaN();
    double TmaxK = std::numeric_limits<double>::quiet_NaN();

    bool hasData = false;
};

class ChemsepComponent {
public:
    void setField(
        const std::string& name,
        const std::string& value
        );

    void setScalar(
        const std::string& name,
        double value
        );

    void setCorrelation(
        const std::string& name,
        const ChemsepCorrelation& correlation
        );

    const std::string& field(
        const std::string& name
        ) const;

    bool hasField(
        const std::string& name
        ) const;

    bool hasScalar(
        const std::string& name
        ) const;

    double scalar(
        const std::string& name
        ) const;

    double scalarOrNaN(
        const std::string& name
        ) const;

    bool hasCorrelation(
        const std::string& name
        ) const;

    const ChemsepCorrelation& correlation(
        const std::string& name
        ) const;

    const std::unordered_map<std::string, std::string>& rawFields() const;
    const std::unordered_map<std::string, double>& scalarFields() const;
    const std::unordered_map<std::string, ChemsepCorrelation>& correlations() const;

    std::string componentName() const;
    std::string formulaChemsep() const;
    std::string casNumber() const;

    double molarMassKgPerKmol() const;
    double criticalTemperatureK() const;
    double criticalPressurePa() const;
    double criticalVolumeM3PerKmol() const;
    double criticalCompressibility() const;
    double acentricFactor() const;
    double idealGasFormationEnthalpyJPerKmol() const;

private:
    std::unordered_map<std::string, std::string> rawFields_;
    std::unordered_map<std::string, double> scalarFields_;
    std::unordered_map<std::string, ChemsepCorrelation> correlations_;
};
