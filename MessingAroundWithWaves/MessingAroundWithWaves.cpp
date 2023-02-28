#include "ExplicitWorld2d.h"
#include "Display.h"

/* Sources:
https://www.simscale.com/docs/simwiki/numerics-background/what-are-boundary-conditions/
https://youtu.be/MqlDQzw3x6Q
https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
https://www.youtube.com/watch?v=6L4Ok_zvZ6c
*/

// Explicit world parameters.
static constexpr const double FPS = 240.0f;
static constexpr const double RESOLUTION = 20.0f;
static constexpr const double CELERITY = 343.0f * 0.01f; // ~1/3 of a meter per second
static constexpr const double CFL = 2.0f * (CELERITY * (1.0f / FPS) / (1.0f / RESOLUTION));

// Source parameters.
static constexpr const double Pi = 3.14159265359f;
static constexpr const double FREQUENCY = 6.384f * 1.f; // 3 hertz
static constexpr const double AMPLITUDE = 1.0f;
double GenerateSourceSample(const double time) {
    return AMPLITUDE * std::cos(2.0f * Pi * time * FREQUENCY);
}

// Display parameters.
static constexpr const size_t DISPLAY_RESOLUTION = 600;

int main(void)
{
    DisplacementField2D initialConditions((size_t)RESOLUTION, std::vector<double>((size_t)RESOLUTION, 0.0f)); // Undisturbed initial conditions.
    std::vector<Source> sources = { Source(GenerateSourceSample, {0.5f, 0.75f}) }; // Source at (0.5, 0.75)
    std::vector<Rect> rects = { Rect({0.25f, 0.25f},{0.75f, 0.5f}) }; // Rect from (0.25, 0.25) to (0.75, 0.5)

    ExplicitWorld2D w(CELERITY, initialConditions, FieldUpdate::HAROON, BoundaryCondition::MUR, sources, rects);

    InitWindow((int)DISPLAY_RESOLUTION, (int)DISPLAY_RESOLUTION, nullptr);
    SetTargetFPS((int)FPS);

    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        BeginDrawing();
        ClearBackground(RAYWHITE);

        DrawDisplacement(w.Update(1.0 / FPS), DISPLAY_RESOLUTION, DISPLAY_RESOLUTION);

        EndDrawing();
    }

    CloseWindow();        // Close window and OpenGL context
    return 0;
}