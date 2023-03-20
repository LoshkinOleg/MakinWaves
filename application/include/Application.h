#pragma once

#include <cmath>

#include "raylib.h"

#include "ExplicitWorld2d.h"

namespace Sim {

    /*
        Utility function to remap a value from one range to another. Used to convert the [-inf ; +inf] range in which the displacement can be to a [0 ; 1] range for display with raylib.
    */
    constexpr inline const float RemapToRange(const float inRangeMin, const float inRangeMax, const float outRangeMin, const float outRangeMax, const float value)
    {
        return outRangeMin + (value - inRangeMin) * (outRangeMax - outRangeMin) / (inRangeMax - inRangeMin);
    }

    /*
        Utility function to help us see the tiny differences in pression.
        Desmos visualization of the function: https://www.desmos.com/calculator/6n7u2lzoyu
    */
    inline float ContrastSigmoid(const float value, const float contrast)
    {
        return (1.0f / (1.0f + std::expf(-value * contrast))) - 0.5f;
    }
    constexpr const float CONTRAST = 1.0f;

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
                        const auto contrastedDisplacement = ContrastSigmoid(displacement, CONTRAST);
                        const auto remappedDisplacement = RemapToRange(-1.0, 1.0, 0.0, 1.0, contrastedDisplacement);
                        const unsigned char displacementForRaylib = std::clamp(remappedDisplacement, 0.0f, 1.0f) * 255;
                        DrawPixel(displayX, displayY,
                            { displacementForRaylib , displacementForRaylib , displacementForRaylib , 255 } // {r,g,b,a} with unsigned chars.
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

        ExplicitWorld2d w_;
        DisplacementField2D dToDraw_;
        const size_t simResX_, simResY_;
    };
}