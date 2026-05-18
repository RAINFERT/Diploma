#include "NewtonSolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace numerics
{
void NewtonSolver::validateVector(
    const Vector& vector,
    const std::string& name
    )
{
    if (vector.empty())
    {
        throw std::invalid_argument(name + " must not be empty");
    }

    for (double value : vector)
    {
        if (!std::isfinite(value))
        {
            throw std::invalid_argument(name + " contains non-finite value");
        }
    }
}

void NewtonSolver::validateOptions(
    const NewtonOptions& options,
    std::size_t dimension
    )
{
    if (options.maxIterations <= 0)
    {
        throw std::invalid_argument("Newton maxIterations must be positive");
    }

    if (options.residualTolerance <= 0.0)
    {
        throw std::invalid_argument("Newton residualTolerance must be positive");
    }

    if (options.stepTolerance <= 0.0)
    {
        throw std::invalid_argument("Newton stepTolerance must be positive");
    }

    if (options.finiteDifferenceStep <= 0.0)
    {
        throw std::invalid_argument("Newton finiteDifferenceStep must be positive");
    }

    if (options.minDampingFactor <= 0.0 || options.minDampingFactor > 1.0)
    {
        throw std::invalid_argument("Newton minDampingFactor must be in (0, 1]");
    }

    if (!options.lowerBounds.empty() && options.lowerBounds.size() != dimension)
    {
        throw std::invalid_argument("Newton lowerBounds size mismatch");
    }

    if (!options.upperBounds.empty() && options.upperBounds.size() != dimension)
    {
        throw std::invalid_argument("Newton upperBounds size mismatch");
    }

    if (!options.lowerBounds.empty() && !options.upperBounds.empty())
    {
        for (std::size_t i = 0; i < dimension; ++i)
        {
            if (options.lowerBounds[i] > options.upperBounds[i])
            {
                throw std::invalid_argument("Newton lower bound is greater than upper bound");
            }
        }
    }
}

bool NewtonSolver::isInsideBounds(
    const Vector& x,
    const NewtonOptions& options
    )
{
    if (!options.lowerBounds.empty())
    {
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            if (x[i] < options.lowerBounds[i])
            {
                return false;
            }
        }
    }

    if (!options.upperBounds.empty())
    {
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            if (x[i] > options.upperBounds[i])
            {
                return false;
            }
        }
    }

    return true;
}

double NewtonSolver::infinityNorm(
    const Vector& vector
    )
{
    double result = 0.0;

    for (double value : vector)
    {
        result = std::max(result, std::abs(value));
    }

    return result;
}

Vector NewtonSolver::addScaled(
    const Vector& x,
    const Vector& dx,
    double scale
    )
{
    if (x.size() != dx.size())
    {
        throw std::invalid_argument("Vector size mismatch in addScaled");
    }

    Vector result(x.size());

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        result[i] = x[i] + scale * dx[i];
    }

    return result;
}

Matrix NewtonSolver::finiteDifferenceJacobian(
    const Function& function,
    const Vector& x,
    const Vector& functionValue,
    double step
    )
{
    validateVector(x, "Jacobian x");
    validateVector(functionValue, "Jacobian functionValue");

    if (step <= 0.0)
    {
        throw std::invalid_argument("Finite difference step must be positive");
    }

    const std::size_t n = x.size();

    if (functionValue.size() != n)
    {
        throw std::invalid_argument("Newton requires square system: F dimension must match x dimension");
    }

    Matrix jacobian(
        n,
        Vector(n, 0.0)
        );

    for (std::size_t column = 0; column < n; ++column)
    {
        Vector xPerturbed = x;

        const double perturbation =
            step * std::max(1.0, std::abs(x[column]));

        xPerturbed[column] += perturbation;

        Vector fPerturbed =
            function(xPerturbed);

        if (fPerturbed.size() != n)
        {
            throw std::runtime_error("Function returned vector with invalid size during Jacobian calculation");
        }

        validateVector(fPerturbed, "Jacobian fPerturbed");

        for (std::size_t row = 0; row < n; ++row)
        {
            jacobian[row][column] =
                (fPerturbed[row] - functionValue[row])
                / perturbation;
        }
    }

    return jacobian;
}

Vector NewtonSolver::solveLinearSystem(
    Matrix matrix,
    Vector rightHandSide
    )
{
    const std::size_t n = matrix.size();

    if (n == 0)
    {
        throw std::invalid_argument("Linear system matrix must not be empty");
    }

    if (rightHandSide.size() != n)
    {
        throw std::invalid_argument("Linear system right-hand side size mismatch");
    }

    for (const Vector& row : matrix)
    {
        if (row.size() != n)
        {
            throw std::invalid_argument("Linear system matrix must be square");
        }
    }

    constexpr double PivotTolerance = 1.0e-14;

    // Прямой ход метода Гаусса с частичным выбором главного элемента.
    for (std::size_t column = 0; column < n; ++column)
    {
        std::size_t pivotRow = column;
        double pivotAbs = std::abs(matrix[column][column]);

        for (std::size_t row = column + 1; row < n; ++row)
        {
            const double candidateAbs =
                std::abs(matrix[row][column]);

            if (candidateAbs > pivotAbs)
            {
                pivotAbs = candidateAbs;
                pivotRow = row;
            }
        }

        if (pivotAbs < PivotTolerance)
        {
            throw std::runtime_error("Linear system is singular or ill-conditioned");
        }

        if (pivotRow != column)
        {
            std::swap(matrix[pivotRow], matrix[column]);
            std::swap(rightHandSide[pivotRow], rightHandSide[column]);
        }

        const double pivot =
            matrix[column][column];

        for (std::size_t row = column + 1; row < n; ++row)
        {
            const double factor =
                matrix[row][column] / pivot;

            matrix[row][column] = 0.0;

            for (std::size_t col = column + 1; col < n; ++col)
            {
                matrix[row][col] -=
                    factor * matrix[column][col];
            }

            rightHandSide[row] -=
                factor * rightHandSide[column];
        }
    }

    // Обратный ход.
    Vector solution(n, 0.0);

    for (std::size_t rowReverse = 0; rowReverse < n; ++rowReverse)
    {
        const std::size_t row =
            n - 1 - rowReverse;

        double sum =
            rightHandSide[row];

        for (std::size_t col = row + 1; col < n; ++col)
        {
            sum -=
                matrix[row][col] * solution[col];
        }

        const double diagonal =
            matrix[row][row];

        if (std::abs(diagonal) < PivotTolerance)
        {
            throw std::runtime_error("Linear system has zero diagonal during back substitution");
        }

        solution[row] =
            sum / diagonal;
    }

    return solution;
}

NewtonResult NewtonSolver::solve(
    const Function& function,
    const Vector& initialGuess,
    const NewtonOptions& options
    )
{
    validateVector(initialGuess, "Newton initialGuess");
    validateOptions(options, initialGuess.size());

    if (!isInsideBounds(initialGuess, options))
    {
        throw std::invalid_argument("Newton initialGuess is outside bounds");
    }

    NewtonResult result;
    result.solution = initialGuess;

    Vector x = initialGuess;

    for (int iteration = 0; iteration <= options.maxIterations; ++iteration)
    {
        Vector fx =
            function(x);

        if (fx.size() != x.size())
        {
            throw std::runtime_error("Newton function returned vector with invalid size");
        }

        validateVector(fx, "Newton residual");

        const double residualNorm =
            infinityNorm(fx);

        result.iterations = iteration;
        result.solution = x;
        result.residual = fx;
        result.residualNorm = residualNorm;

        if (residualNorm < options.residualTolerance)
        {
            result.converged = true;
            result.message = "Converged by residual norm";
            return result;
        }

        if (iteration == options.maxIterations)
        {
            break;
        }

        Matrix jacobian =
            finiteDifferenceJacobian(
                function,
                x,
                fx,
                options.finiteDifferenceStep
                );

        Vector minusFx(fx.size());

        for (std::size_t i = 0; i < fx.size(); ++i)
        {
            minusFx[i] = -fx[i];
        }

        Vector dx =
            solveLinearSystem(
                jacobian,
                minusFx
                );

        const double fullStepNorm =
            infinityNorm(dx);

        result.stepNorm = fullStepNorm;

        if (fullStepNorm < options.stepTolerance)
        {
            result.converged = true;
            result.message = "Converged by step norm";
            return result;
        }

        double damping = 1.0;
        bool accepted = false;

        Vector bestCandidate = x;
        Vector bestResidual = fx;
        double bestResidualNorm = residualNorm;

        while (damping >= options.minDampingFactor)
        {
            Vector candidate =
                addScaled(
                    x,
                    dx,
                    damping
                    );

            if (!isInsideBounds(candidate, options))
            {
                damping *= 0.5;
                continue;
            }

            Vector candidateResidual =
                function(candidate);

            if (candidateResidual.size() != x.size())
            {
                throw std::runtime_error("Newton function returned invalid vector size during line search");
            }

            validateVector(candidateResidual, "Newton candidateResidual");

            const double candidateResidualNorm =
                infinityNorm(candidateResidual);

            if (candidateResidualNorm < bestResidualNorm)
            {
                bestCandidate = candidate;
                bestResidual = candidateResidual;
                bestResidualNorm = candidateResidualNorm;
            }

            // Принимаем шаг, если он уменьшил норму остатка.
            if (candidateResidualNorm < residualNorm)
            {
                x = candidate;
                accepted = true;
                break;
            }

            damping *= 0.5;
        }

        if (!accepted)
        {
            // Если ни один демпфированный шаг не прошел строгий критерий,
            // но лучший найденный вариант все-таки улучшает residual,
            // используем его.
            if (bestResidualNorm < residualNorm)
            {
                x = bestCandidate;
                accepted = true;
            }
        }

        if (!accepted)
        {
            result.converged = false;
            result.solution = x;
            result.residual = fx;
            result.residualNorm = residualNorm;
            result.message = "Line search failed to reduce residual";
            return result;
        }
    }

    result.converged = false;
    result.message = "Maximum iterations reached";
    return result;
}
}
