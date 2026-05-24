#pragma once

#include <cstddef>
#include <stdexcept>

class DynamicStateLayout {
public:
    explicit DynamicStateLayout(std::size_t componentCount)
        : componentCount_(componentCount)
    {
        if (componentCount_ == 0) {
            throw std::runtime_error(
                "DynamicStateLayout: component count must be positive"
                );
        }
    }

    std::size_t componentCount() const {
        return componentCount_;
    }

    // ============================================================
    // Variable indices
    // ============================================================

    std::size_t mass(std::size_t componentIndex) const {
        checkComponentIndex(componentIndex);
        return componentIndex;
    }

    std::size_t energy() const {
        return componentCount_;
    }

    std::size_t pressure() const {
        return componentCount_ + 1;
    }

    std::size_t temperature() const {
        return componentCount_ + 2;
    }

    std::size_t liquidSplit(std::size_t componentIndex) const {
        checkComponentIndex(componentIndex);
        return componentCount_ + 3 + componentIndex;
    }

    std::size_t liquidMolarVolume() const {
        return 2 * componentCount_ + 3;
    }

    std::size_t vaporMolarVolume() const {
        return 2 * componentCount_ + 4;
    }

    std::size_t variableCount() const {
        return 2 * componentCount_ + 5;
    }

    // ============================================================
    // Residual indices
    // ============================================================

    std::size_t massResidual(std::size_t componentIndex) const {
        checkComponentIndex(componentIndex);
        return componentIndex;
    }

    std::size_t energyResidual() const {
        return componentCount_;
    }

    std::size_t volumeResidual() const {
        return componentCount_ + 1;
    }

    std::size_t inventoryEnergyResidual() const {
        return componentCount_ + 2;
    }

    std::size_t liquidSplitOdeResidual(std::size_t componentIndex) const {
        checkComponentIndex(componentIndex);
        return componentCount_ + 3 + componentIndex;
    }

    std::size_t liquidMolarVolumeOdeResidual() const {
        return 2 * componentCount_ + 3;
    }

    std::size_t vaporMolarVolumeOdeResidual() const {
        return 2 * componentCount_ + 4;
    }

    std::size_t residualCount() const {
        return 2 * componentCount_ + 5;
    }

private:
    std::size_t componentCount_;

    void checkComponentIndex(std::size_t componentIndex) const {
        if (componentIndex >= componentCount_) {
            throw std::runtime_error(
                "DynamicStateLayout: component index out of range"
                );
        }
    }
};
