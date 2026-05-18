#pragma once

#include <functional>
#include <string>
#include <vector>

namespace numerics
{
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

struct NewtonOptions
{
    int maxIterations = 50;

    double residualTolerance = 1.0e-8;
    double stepTolerance = 1.0e-10;

    // Шаг для численного расчета Якоби
    double finiteDifferenceStep = 1.0e-6;

    // Минимальный коэффициент демпфирования
    double minDampingFactor = 1.0e-6;

    // Если векторы пустые, ограничения не используются.
    Vector lowerBounds;
    Vector upperBounds;
};

struct NewtonResult
{
    bool converged = false;

    int iterations = 0;

    Vector solution;
    Vector residual;

    double residualNorm = 0.0;
    double stepNorm = 0.0;

    std::string message;
};

class NewtonSolver
{
public:
    using Function = std::function<Vector(const Vector&)>;

    static NewtonResult solve(
        const Function& function,
        const Vector& initialGuess,
        const NewtonOptions& options = NewtonOptions()
        );

    static Matrix finiteDifferenceJacobian(
        const Function& function,
        const Vector& x,
        const Vector& functionValue,
        double step
        );

    static Vector solveLinearSystem(
        Matrix matrix,
        Vector rightHandSide
        );

    static double infinityNorm(
        const Vector& vector
        );

private:
    static void validateVector(
        const Vector& vector,
        const std::string& name
        );

    static void validateOptions(
        const NewtonOptions& options,
        std::size_t dimension
        );

    static bool isInsideBounds(
        const Vector& x,
        const NewtonOptions& options
        );

    static Vector addScaled(
        const Vector& x,
        const Vector& dx,
        double scale
        );
};
}
