#include <stdexcept>
#include <vector>
#include <functional>
#include <algorithm>

#include <raylib.h>

namespace Sim {
    constexpr const float Pi = 3.14159265359;

    /*
        // An undisturbed displacement one that is equal to 0.0f .
        TODO: decribe
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
    
    class Source {
    public:
        Source() = delete;
        Source(const std::function<float(float)> sampleProvidingCallback, const float normalizedPosX, const float normalizedPosY):
            callback_(sampleProvidingCallback), posX(normalizedPosX), posY(normalizedPosY) {}

        inline float GetSample(const float t) const {
            return callback_(t);
        }

        float posX;
        float posY;

    private:
        const std::function<float(float)> callback_;
    };

    struct Obstacle {
    public:
        Obstacle(const float normalizedMinX, const float normalizedMinY, const float normalizedMaxX, const float normalizedMaxY) :
            minX(normalizedMinX), minY(normalizedMaxY), maxX(normalizedMaxX), maxY(normalizedMaxY) {}

        float minX = 0.0f;
        float minY = 0.0f;
        float maxX = 0.0f;
        float maxY = 0.0f;
    };

    /*
        A description of the CFL: https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
        TODO: explain
    */
    constexpr inline float CourantFreidrichsLevyCondition(const float celerityX, const float celerityY, const float dt, const float dx, const float dy) {
        return (celerityX * dt / dx) + (celerityY * dt / dy);
    }

    constexpr inline const float RemapToRange(const float inRangeMin, const float inRangeMax, const float outRangeMin, const float outRangeMax, const float value)
    {
        return outRangeMin + (value - inRangeMin) * (outRangeMax - outRangeMin) / (inRangeMax - inRangeMin);
    }

    inline float ContrastSigmoid(const float x, const float f)
    {
        return (1.0f / (1.0f + std::expf(-x * f))) - 0.5f;
    }

    /*
        TODO: description + link to desmos graph
    */
    inline float GenerateSourceSample(const float t, const float frq, const float amp) {
        return amp * std::cos(2.0f * Pi * t * frq);
    }

    enum BoundaryCondition {
        DIRICHLET = 0, // Displacement at boundaries is always 0. Causes perfect reflections (no loss of energy).
        MUR = 1 // Mur's boundary condition. Absorbs most of the energy (but not all).
    };

    class ExplicitWorld2d {
    public:
        ExplicitWorld2d() = delete;
        ExplicitWorld2d(
            const float deltaTime,
            const float waveCelerityX, const float waveCelerityY,
            const size_t simulationResolutionX, const size_t simulationResolutionY,
            const BoundaryCondition boundaryConditionToUse,
            const std::vector<Source>& sources, const std::vector<Obstacle>& obstacles)
            :
            dLast_(DisplacementField2D(simulationResolutionX, simulationResolutionY)), dCurrent_(DisplacementField2D(simulationResolutionX, simulationResolutionY)), dNext_(DisplacementField2D(simulationResolutionX, simulationResolutionY)),
            sources_(sources), obstacles_(obstacles),
            cfl_(CourantFreidrichsLevyCondition(waveCelerityX, waveCelerityY, deltaTime, 1.0f / simulationResolutionX, 1.0f / simulationResolutionY)),
            simResX_(simulationResolutionX),
            simResY_(simulationResolutionY),
            boundaryCond_(boundaryConditionToUse),
            dt_(deltaTime)
        {
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

                for (size_t x = 1; x < simResX_ - 1; x++) // Edges of the field are boundaries, so they will be updated using specified boundary conditions, hence the -1 and 1 .
                {
                    for (size_t y = 1; y < simResY_ - 1; y++)
                    {
                        const float currentCenter =     dCurrent_[x + 0][y + 0];
                        const float currentRight =      dCurrent_[x + 1][y + 0];
                        const float currentLeft =       dCurrent_[x - 1][y + 0];
                        const float currentTop =        dCurrent_[x + 0][y - 1];
                        const float currentBottom =     dCurrent_[x + 0][y + 1];

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

                t_++;
            }

            std::copy(dCurrent_.d.begin(), dCurrent_.d.end(), dLast_.d.begin());
            std::copy(dNext_.d.begin(), dNext_.d.end(), dCurrent_.d.begin());
            return dNext_;
        }

    private:
        DisplacementField2D dLast_, dCurrent_, dNext_;
        std::vector<Source> sources_;
        std::vector<Obstacle> obstacles_;
        size_t t_ = 0.0f;
        float accumulatedDt_ = 0.0f;
        const float cfl_;
        const size_t simResX_, simResY_;
        const BoundaryCondition boundaryCond_;
        const float dt_;
    };

    class Application {
    public:
        Application() = delete;
        Application(
            const size_t targetDisplayFPS,
            const size_t displayResolutionX, const size_t displayResolutionY,
            const float deltaTime,
            const float waveCelerityX, const float waveCelerityY,
            const size_t simulationResolutionX, const size_t simulationResolutionY,
            const BoundaryCondition boundaryConditionToUse,
            const std::vector<Source>& sources, const std::vector<Obstacle>& obstacles)
            :
            DisplayRefreshRate(targetDisplayFPS),
            DisplayResolutionX(displayResolutionX), DisplayResolutionY(displayResolutionY),
            simResX_(simulationResolutionX), simResY_(simulationResolutionY),
            dToDraw_(DisplacementField2D(simulationResolutionX, simulationResolutionY)),
            w_(ExplicitWorld2d(deltaTime, waveCelerityX, waveCelerityY, simulationResolutionX, simulationResolutionY, boundaryConditionToUse, sources, obstacles))
        {}

        void Run(const float timeScale) {
            InitWindow((int)DisplayResolutionX, (int)DisplayResolutionY, nullptr);
            SetTargetFPS((int)DisplayRefreshRate);

            while (!WindowShouldClose())    // Detect window close button or ESC key
            {
                BeginDrawing();
                ClearBackground(RAYWHITE);

                for (size_t displayX = 0; displayX < DisplayResolutionX; displayX++)
                {
                    for (size_t displayY = 0; displayY < DisplayResolutionY; displayY++)
                    {
                        const float normalizedCoordX = (float)displayX / (float)DisplayResolutionX; // Normalize our coordinate so that we can match it to the resolution of the simulation.
                        const float normalizedCoordY = (float)displayY / (float)DisplayResolutionY;
                        const size_t simCoordX = std::floor((float)simResX_ * normalizedCoordX); // Convert normalized coordinate to the simulation coordinate.
                        const size_t simCoordY = std::floor((float)simResY_ * normalizedCoordY);
                        const auto displacement = dToDraw_[simCoordX][simCoordY]; // Retrieve the value of the displacement.
                        const auto contrastedDisplacement = ContrastSigmoid(displacement, 1);
                        const auto remappedDisplacement = RemapToRange(-1.0, 1.0, 0.0, 1.0, contrastedDisplacement);
                        const unsigned char displacementForRaylib = std::clamp(remappedDisplacement, 0.0f, 1.0f) * 255;
                        DrawPixel(displayX, displayY,
                            { displacementForRaylib , displacementForRaylib , displacementForRaylib , 255} // {r,g,b,a} with unsigned chars.
                        );
                    }
                }

                EndDrawing();

                const auto& dNext = w_.Update(GetFrameTime() * timeScale);
                std::copy(dNext.d.begin(), dNext.d.end(), dToDraw_.d.begin());
            }

            CloseWindow();
        }

        const size_t DisplayRefreshRate;
        const size_t DisplayResolutionX, DisplayResolutionY;
    private:

        const size_t simResX_, simResY_;
        DisplacementField2D dToDraw_;
        ExplicitWorld2d w_;
    };
}

/* Sources:
https://www.simscale.com/docs/simwiki/numerics-background/what-are-boundary-conditions/
https://youtu.be/MqlDQzw3x6Q
https://www.youtube.com/watch?v=6L4Ok_zvZ6c
https://www.youtube.com/watch?v=6xwLRfYxCno
*/

int main(void)
{
    const auto sourceCallback = [](float t)->float {return Sim::GenerateSourceSample(t, 20, 1000000000); };

    Sim::Source source = Sim::Source(sourceCallback, 0.5f, 0.5f);
    Sim::Application app = Sim::Application(30, 600, 600, 0.0001, 24, 24, 100, 100, Sim::MUR, {source}, {});
    app.Run(1000);
    return 0;
}