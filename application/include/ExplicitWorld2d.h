#pragma once

#include <stdexcept>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <string>

namespace Sim {
    /*
        The Courant–Friedrichs–Lewy condition, a formula used in the computations below and one that allows you to detect when a wave moves across the simulation grid faster than 1 cell / update, which would break the simulation.
        A description of the CFL: https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
    */
    constexpr inline float CourantFreidrichsLevyCondition(const float celerity, const float dt, const float dx) {
        return 2.0f * (celerity * dt / dx);
    }

    enum BoundaryCondition {
        DIRICHLET = 0, // Displacement at boundaries is always 0. Causes perfect reflections (no loss of energy).
        MUR = 1 // Mur's boundary condition. Absorbs most of the energy (but not all).
    };

    /*
        Just a 2D array of floats to represent the displacement at a particular position.
        In our case, an undisturbed displacement is one that is equal to 0.0f.
    */
    struct DisplacementField2D
    {
        DisplacementField2D(const size_t spatialResolutionX, const size_t spatialResolutionY) {
            d = std::vector<std::vector<float>>(spatialResolutionX, std::vector<float>(spatialResolutionY, 0.0f));
        }

        inline std::vector<float>& operator[](const size_t idx) {
            return d[idx];
        }

        std::vector<std::vector<float>> d;
    };

    /*
        Just a class that will provide the ExlpicitWorld2d with a sample.
        Has a normalized position and that's about it.
    */
    class Source {
    public:
        Source() = delete;
        Source(const std::function<float(float)> sampleProvidingCallback, const float normalizedPosX, const float normalizedPosY) :
            callback_(sampleProvidingCallback), posX(normalizedPosX), posY(normalizedPosY) {}

        inline float GetSample(const float t) const {
            return callback_(t);
        }

        float posX;
        float posY;

    private:
        const std::function<float(float)> callback_;
    };

    /*
        Just a struct that will provide the ExplicitWorld2d with the position of an obstacle.
        Has a normalized min and max coordinates and that's it.
    */
    struct Obstacle {
    public:
        Obstacle(const float normalizedMinX, const float normalizedMinY, const float normalizedMaxX, const float normalizedMaxY) :
            minX(normalizedMinX), minY(normalizedMinY), maxX(normalizedMaxX), maxY(normalizedMaxY) {}

        float minX = 0.0f;
        float minY = 0.0f;
        float maxX = 0.0f;
        float maxY = 0.0f;
    };

    /*
        The meat and gravy of the simulation code. Updates and tracks the last, current and the next displacement of the world.
        It's called explicit because if we want to compute the state of the world at a time T, we'd first need to compute the state of the world at all the previous steps, in order. There exist implicit techniques for solving the same issue where that's not a concern.
    */
    class ExplicitWorld2d {
    public:
        ExplicitWorld2d() = delete;
        ExplicitWorld2d(
            const float deltaTime,
            const float waveCelerity,
            const size_t simulationResolution,
            const BoundaryCondition boundaryConditionToUse,
            const std::vector<Source>& sources, const std::vector<Obstacle>& obstacles)
            :
            dLast_(DisplacementField2D(simulationResolution, simulationResolution)), dCurrent_(DisplacementField2D(simulationResolution, simulationResolution)), dNext_(DisplacementField2D(simulationResolution, simulationResolution)),
            sources_(sources), obstacles_(obstacles),
            cfl_(CourantFreidrichsLevyCondition(waveCelerity, deltaTime, 1.0f / simulationResolution)),
            simResX_(simulationResolution),
            simResY_(simulationResolution),
            boundaryCond_(boundaryConditionToUse),
            dt_(deltaTime)
        {
            std::cout << "CFL value: " << std::to_string(cfl_) << std::endl;
            if (cfl_ > 1.0f) {
                throw std::runtime_error("CourantFreidrichsLevy condition is not fulfilled, meaning that the wave propagates too quickly for the chosen grid size and the simulation is therefore not accurate.");
            }
        }

        const DisplacementField2D& Update(const float timeSinceLastDisplayFrame)
        {
            // Update only if enough real time has passed to warrant a simulation update, and update multiple times if more than one dt time has passed since the last display update.
            accumulatedDt_ += timeSinceLastDisplayFrame;
            const size_t iterations = std::floor(accumulatedDt_);
            accumulatedDt_ -= iterations;
            for (size_t it = 0; it < iterations; it++)
            {
                /*
                Neighboring cells:
                       |  X - 1  |   X     | X + 1
                 ----- | --------------------------
                 Y - 1 | -       | Top     | -
                     Y | Left    | Center  | Right
                 Y + 1 | -       | Bottom  | -
                */
                // dNextCenter = 2*dCurrentCenter - dLastCenter + CFL^2 * (dCurrentRight + dCurrentLeft - 4*dCurrentCenter + dCurrentBottom + dCurrentTop) + dt^2 * f(t)

                for (size_t x = 1; x < simResX_ - 1; x++) // Edges of the field are boundaries, so they will be updated using specified boundary conditions, hence the -1 and 1 .
                {
                    for (size_t y = 1; y < simResY_ - 1; y++)
                    {
                        const float currentCenter = dCurrent_[x + 0][y + 0];
                        const float currentRight = dCurrent_[x + 1][y + 0];
                        const float currentLeft = dCurrent_[x - 1][y + 0];
                        const float currentTop = dCurrent_[x + 0][y - 1];
                        const float currentBottom = dCurrent_[x + 0][y + 1];

                        const float lastCenter = dLast_[x][y];

                        // Formula when no source is present (when f(x) is a constant 0):
                        // dNextCenter = 2*dCurrentCenter - dLastCenter + CFL^2 * (dCurrentRight + dCurrentLeft - 4*dCurrentCenter + dCurrentBottom + dCurrentTop)
                        dNext_[x][y] =
                            2.0f * currentCenter
                            - lastCenter
                            + (cfl_ * cfl_) // dt, dx, dy, cx and cy are hidden in there.
                            * (currentRight + currentLeft - 4.0f * currentCenter + currentTop + currentBottom);
                    }
                }

                // Update field with forcing functions (aka sound sources). Doing it in a separate pass to avoid checking for every cell whether it contains a source.
                for (auto& src : sources_)
                {
                    const size_t simPosX = std::clamp(src.posX, 0.0f, 1.0f) * simResX_;
                    const size_t simPosY = std::clamp(src.posY, 0.0f, 1.0f) * simResY_;

                    // Formula: dNextCenter += dt^2 * f(t). This is the part we've left out in the pass above and assumed to be a constant function f(x) = 0 .
                    dNext_[simPosX][simPosY] += dt_ * dt_ * src.GetSample(t_ * dt_);
                }

                // Finally, update boundary conditions.
                switch (boundaryCond_)
                {
                case BoundaryCondition::DIRICHLET: // Essentially, we're saying that displacement at the boundaries should always be at 0 .
                {
                    // Left boundary. x = left boundary, y varies.
                    for (size_t y = 0; y < simResY_; y++)
                    {
                        // Formula: dNextCenter = 0
                        dNext_[0][y] = 0.0f;
                    }

                    // Right boundary. x = right boundary, y varies.
                    for (size_t y = 0; y < simResY_; y++)
                    {
                        // Formula: dNextCenter = 0
                        dNext_[simResX_ - 1][y] = 0.0f;
                    }

                    // Top boundary. x varies, y = top boundary.
                    for (size_t x = 0; x < simResX_; x++)
                    {
                        // Formula: dNextCenter = 0
                        dNext_[x][0] = 0.0f;
                    }

                    // Bottom boundary. x varies, y = bottom boundary.
                    for (size_t x = 0; x < simResX_; x++)
                    {
                        // Formula: dNextCenter = 0
                        dNext_[x][simResY_ - 1] = 0.0f;
                    }

                    // Obstacles.
                    for (auto& obj : obstacles_)
                    {
                        const size_t minPosX = std::clamp(obj.minX, 0.0f, 1.0f) * simResX_;
                        const size_t minPosY = std::clamp(obj.minY, 0.0f, 1.0f) * simResY_;
                        const size_t maxPosX = std::clamp(obj.maxX, 0.0f, 1.0f) * simResX_;
                        const size_t maxPosY = std::clamp(obj.maxY, 0.0f, 1.0f) * simResY_;

                        for (size_t x = minPosX; x < maxPosX; x++)
                        {
                            for (size_t y = minPosY; y < maxPosY; y++)
                            {
                                // Formula: dNextCenter = 0
                                dNext_[x][y] = 0.0f;
                            }
                        }
                    }
                } break;
                case BoundaryCondition::MUR:
                {
                    // Left boundary. x = left boundary, y varies.
                    for (size_t y = 0; y < simResY_; y++)
                    {
                        const double currentRight = dCurrent_[1][y];
                        const double currentCenter = dCurrent_[0][y];
                        const double nextRight = dNext_[1][y];

                        // Formula: dNextCenter = dCurrentRight + (CFL - 1 / CFL + 1) * (dNextRight - dCurrentCenter)
                        dNext_[0][y] =
                            currentRight
                            + ((cfl_ - 1.0f) / (cfl_ + 1.0f))
                            * (nextRight - currentCenter);
                    }
                    // Right boundary. x = right boundary, y varies.
                    for (size_t y = 0; y < simResY_; y++)
                    {
                        const double currentLeft = dCurrent_[simResX_ - 2][y];
                        const double currentCenter = dCurrent_[simResX_ - 1][y];
                        const double nextLeft = dNext_[simResX_ - 2][y];

                        // Formula: dNextCenter = dCurrentLeft + (CFL - 1 / CFL + 1) * (dNextLeft - dCurrentCenter)
                        dNext_[simResX_ - 1][y] =
                            currentLeft
                            + ((cfl_ - 1.0f) / (cfl_ + 1.0f))
                            * (nextLeft - currentCenter);
                    }
                    // Top boundary. x varies, y = top boundary.
                    for (size_t x = 0; x < simResX_; x++)
                    {
                        const double currentBottom = dCurrent_[x][1];
                        const double currentCenter = dCurrent_[x][0];
                        const double nextBottom = dNext_[x][1];

                        // Formula: dNextCenter = dCurrentBottom + (CFL - 1 / CFL + 1) * (dNextBottom - dCurrentCenter)
                        dNext_[x][0] =
                            currentBottom
                            + ((cfl_ - 1.0f) / (cfl_ + 1.0f))
                            * (nextBottom - currentCenter);
                    }
                    // Bottom boundary. x varies, y = bottom boundary.
                    for (size_t x = 0; x < simResX_; x++)
                    {
                        const double currentTop = dCurrent_[x][simResY_ - 2];
                        const double currentCenter = dCurrent_[x][simResY_ - 1];
                        const double nextTop = dNext_[x][simResY_ - 2];

                        // Formula: dNextCenter = dCurrentTop + (CFL - 1 / CFL + 1) * (dNextTop - dCurrentCenter)
                        dNext_[x][simResY_ - 1] =
                            currentTop
                            + ((cfl_ - 1.0f) / (cfl_ + 1.0f))
                            * (nextTop - currentCenter);
                    }

                    // Obstacles.
                    for (auto& obj : obstacles_)
                    {
                        const size_t minPixelPosX = std::clamp(obj.minX, 0.0f, 1.0f) * simResX_;
                        const size_t minPixelPosY = std::clamp(obj.minY, 0.0f, 1.0f) * simResY_;
                        const size_t maxPixelPosX = std::clamp(obj.maxX, 0.0f, 1.0f) * simResX_;
                        const size_t maxPixelPosY = std::clamp(obj.maxY, 0.0f, 1.0f) * simResY_;

                        // Using inverted Mur for the boundaries of the obstacle.

                        // Left boundary. x = left boundary, y varies.
                        for (size_t y = minPixelPosY; y < maxPixelPosY; y++)
                        {
                            const double currentLeft = dCurrent_[minPixelPosX - 1][y];
                            const double currentCenter = dCurrent_[minPixelPosX][y];
                            const double nextLeft = dNext_[minPixelPosX - 1][y];

                            // Formula: dNextCenter = dCurrentLeft + (CFL - 1 / CFL + 1) * (dNextLeft - dCurrentCenter)
                            dNext_[minPixelPosX][y] =
                                currentLeft
                                + ((cfl_ - 1.0f) / (cfl_ + 1.0f))
                                * (nextLeft - currentCenter);
                        }
                        // Right boundary. x = right boundary, y varies.
                        for (size_t y = minPixelPosY; y < maxPixelPosY; y++)
                        {
                            const double currentRight = dCurrent_[maxPixelPosX + 1][y];
                            const double currentCenter = dCurrent_[maxPixelPosX][y];
                            const double nextRight = dNext_[maxPixelPosX + 1][y];

                            // Formula: dNextCenter = dCurrentRight + (CFL - 1 / CFL + 1) * (dNextRight - dCurrentCenter)
                            dNext_[maxPixelPosX][y] =
                                currentRight
                                + ((cfl_ - 1.0f) / (cfl_ + 1.0f))
                                * (nextRight - currentCenter);
                        }
                        // Top boundary. x varies, y = top boundary.
                        for (size_t x = minPixelPosX; x < maxPixelPosX; x++)
                        {
                            const double currentTop = dCurrent_[x][minPixelPosY - 1];
                            const double currentCenter = dCurrent_[x][minPixelPosY];
                            const double nextTop = dNext_[x][minPixelPosY - 1];

                            // Formula: dNextCenter = dCurrentTop + (CFL - 1 / CFL + 1) * (dNextTop - dCurrentCenter)
                            dNext_[x][minPixelPosY] =
                                currentTop
                                + ((cfl_ - 1.0f) / (cfl_ + 1.0f))
                                * (nextTop - currentCenter);
                        }
                        // Bottom boundary. x varies, y = bottom boundary.
                        for (size_t x = minPixelPosX; x < maxPixelPosX; x++)
                        {
                            const double currentBottom = dCurrent_[x][maxPixelPosY + 1];
                            const double currentCenter = dCurrent_[x][maxPixelPosY];
                            const double nextBottom = dNext_[x][maxPixelPosY + 1];

                            // Formula: dNextCenter = dCurrentBottom + (CFL - 1 / CFL + 1) * (dNextBottom - dCurrentCenter)
                            dNext_[x][maxPixelPosY] =
                                currentBottom
                                + ((cfl_ - 1.0f) / (cfl_ + 1.0f))
                                * (nextBottom - currentCenter);
                        }

                        // Using Dirichlet for the unreachable area of the obstacle since the boundaries are supposed to absorb everything.
                        for (size_t x = minPixelPosX + 1; x < maxPixelPosX - 1; x++)
                        {
                            for (size_t y = minPixelPosY + 1; y <= maxPixelPosY - 1; y++)
                            {
                                // Formula: dNextCenter = 0
                                dNext_[x][y] = 0.0f;
                            }
                        }
                    }
                } break;
                default:break;
                }

                std::copy(dCurrent_.d.begin(), dCurrent_.d.end(), dLast_.d.begin());
                std::copy(dNext_.d.begin(), dNext_.d.end(), dCurrent_.d.begin());
                t_++;
            }

            return dNext_;
        }

    private:
        DisplacementField2D dLast_, dCurrent_, dNext_;
        std::vector<Source> sources_;
        std::vector<Obstacle> obstacles_;
        size_t t_ = 0.0f;
        const size_t simResX_, simResY_;
        const BoundaryCondition boundaryCond_;
        float accumulatedDt_ = 0.0f;
        const float cfl_;
        const float dt_;
    };

}