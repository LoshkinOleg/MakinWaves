#pragma once

#include <stdexcept>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <string>

namespace Sim {

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

    // TODO: implement this.
    constexpr inline float CourantFreidrichsLevyCondition(const float celerity, const float dt, const float dx) {
        return 2.0f * celerity * dt / dx;
    }

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

        // TODO: implement this.
        const DisplacementField2D& Update(const float displayDt)
        {
            dtRemainder_ += displayDt;
            const size_t nrOfIterations = std::floor(dtRemainder_);
            dtRemainder_ -= nrOfIterations;

            for (size_t it = 0; it < nrOfIterations; it++)
            {
                /**
                * General formula we're using:
                *   nextCenter =
                *       2*currentAtCenter - lastAtCenter
                *       + (cfl_ * cfl_)
                *       * (currentAtLeft + currentAtTop - 4*currentAtCenter + currentAtRight + currentAtBottom)
                *       + dt_ * dt_ * f(t_ * dt_)
                **/

                for (size_t x = 1; x < simResX_ - 1; x++)
                {
                    for (size_t y = 1; y < simResY_ - 1; y++)
                    {
                        /**
                        *   Neighboring cells:
                        *             |  X - 1  |   X     | X + 1
                        *       ----- | --------------------------
                        *       Y - 1 | -       | Top     | -
                        *           Y | Left    | Center  | Right
                        *       Y + 1 | -       | Bottom  | -
                        **/

                        float& nextCenter = dNext_[x + 0][y + 0];
                        const float currentAtCenter = dCurrent_[x + 0][y + 0];
                        const float lastAtCenter = dLast_[x + 0][y + 0];
                        const float currentAtLeft = dCurrent_[x - 1][y + 0];
                        const float currentAtTop = dCurrent_[x + 0][y - 1];
                        const float currentAtRight = dCurrent_[x + 1][y + 0];
                        const float currentAtBottom = dCurrent_[x + 0][y + 1];

                        nextCenter =
                            2.0f * currentAtCenter - lastAtCenter
                            + (cfl_ * cfl_)
                            * (currentAtLeft + currentAtTop - 4 * currentAtCenter + currentAtRight + currentAtBottom);
                    }
                }

                for (const auto& src : sources_)
                {
                    const size_t simPosX = std::floor(std::clamp(src.posX * (float)simResX_, 1.0f, (float)simResX_ - 1.0f));
                    const size_t simPosY = std::floor(std::clamp(src.posY * (float)simResY_, 1.0f, (float)simResY_ - 1.0f));
                    dNext_[simPosX][simPosY] += dt_ * dt_ * src.GetSample(t_ * dt_);
                }

                switch (boundaryCond_)
                {
                    case DIRICHLET:
                    {
                        for (const auto& obs : obstacles_)
                        {
                            for (size_t x = obs.minX * simResX_; x < obs.maxX * simResX_; x++)
                            {
                                for (size_t y = obs.minY * simResY_; y < obs.maxY * simResY_; y++)
                                {
                                    dNext_[x][y] = 0.0f;
                                }
                            }
                        }

                        for (size_t x = 0; x < simResX_; x++)
                        {
                            const size_t y = 0;
                            dNext_[x][y] = 0.0f;
                        }
                        for (size_t x = 0; x < simResX_; x++)
                        {
                            const size_t y = simResY_ - 1;
                            dNext_[x][y] = 0.0f;
                        }
                        for (size_t y = 0; y < simResY_; y++)
                        {
                            const size_t x = 0;
                            dNext_[x][y] = 0.0f;
                        }
                        for (size_t y = 0; y < simResY_; y++)
                        {
                            const size_t x = simResX_ - 1;
                            dNext_[x][y] = 0.0f;
                        }
                    }break;
                    
                    case MUR:
                    {
                    
                    }break;
                    
                    default:
                        break;
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
        float dtRemainder_ = 0.0f;
        const float cfl_;
        const float dt_;
    };

}