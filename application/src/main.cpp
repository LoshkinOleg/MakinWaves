#include "Application.h"

/* Sources:
https://www.simscale.com/docs/simwiki/numerics-background/what-are-boundary-conditions/
https://youtu.be/MqlDQzw3x6Q
https://www.youtube.com/watch?v=6L4Ok_zvZ6c
https://www.youtube.com/watch?v=6xwLRfYxCno
*/

/*
    A sine function to generate samples for a source. Produces a cosine at the specified frequency and of the specified amplitude.
    Link to graph for visualization: https://www.desmos.com/calculator/az8l1t6bqn 
*/
inline float GenerateSourceSample(const float t, const float frq, const float amp) {
    return amp * std::cos(2.0f * 3.14159265359 * t * frq);
}

// Display parameters.
constexpr const size_t TARGET_DISPLAY_FPS = 30; // Update rate of the display.
constexpr const size_t DISPLAY_RESOLUTION_X = 600; // Resolution of the display.
constexpr const size_t DISPLAY_RESOLUTION_Y = 600;

// Simulation parameters. If you're in Release mode and you get an exception, change to Debug mode and you'll see that it's because you need to adjust the parameters below.
constexpr const float DT = 0.0001f; // Time precision of the simulation. Represents how much time passes with each update of the simulation.
constexpr const float TIME_SCALE = 1000.0f; // Adjust this to make time speed up / slow down depending on your PC's power. The slower the time scale, the easier on the CPU.
constexpr const float WAVE_CELERITY = 1; // Speed at which a wave propagates in the X coordinate (speed of sound). Too slow and the simulation breaks down since the source oscilates faster than the speed of sound. Too fast and the simulation breaks down as well. 24 seems to be the sweet spot for this setup.
constexpr const size_t SIMULATION_GRID_RESOLUTION = 100; // Number of cells on the X coordinate. The bigger the number of cells, the better the resolution, but the smaller DT needs to be prevent the waves from traveling at more than 1 cell per update, which would break the simulation.
constexpr const Sim::BoundaryCondition BOUNDARY_CONDITION = Sim::BoundaryCondition::MUR; // Defines what happens when a wave encouters a physical obstacle. DIRICHLET means that the wave is entirely reflected, MUR means that the wave is (almost) totally absorbed.

// Source parameters. The source is a sine wave.
constexpr const float SOURCE_FREQUENCY = 20.0f; // In hertz.
constexpr const float SOURCE_AMPLITUDE = 1000000.0f; // Needs to be very big to see anything due to the dt^2 factor when computing the effect of the forcing function.
constexpr const float SOURCE_POS_X = 0.5f; // Position of the source on the screen / simulation world as represented with a value ranging [0.0f ; 1.0f]
constexpr const float SOURCE_POS_Y = 0.85f;

// Obstacle parameters. This is a rectangular block.
constexpr const float OBSTACLE_MIN_X = 0.25f; // Top left corner of the rectangular obstacle.
constexpr const float OBSTACLE_MIN_Y = 0.25f;
constexpr const float OBSTACLE_MAX_X = 0.75f; // Bottom right corner of the rectangular obstacle.
constexpr const float OBSTACLE_MAX_Y = 0.50f;

int main(void)
{
    const auto sourceCallback = [](float t)->float {return GenerateSourceSample(t, SOURCE_FREQUENCY, SOURCE_AMPLITUDE); };

    Sim::Source source = Sim::Source(sourceCallback, SOURCE_POS_X, SOURCE_POS_Y);

    Sim::Obstacle obstacle = Sim::Obstacle(OBSTACLE_MIN_X, OBSTACLE_MIN_Y, OBSTACLE_MAX_X, OBSTACLE_MAX_Y);
    
    Sim::Application app =
        Sim::Application(
            TARGET_DISPLAY_FPS, DISPLAY_RESOLUTION_X, DISPLAY_RESOLUTION_Y, // Display parameters.
            DT, WAVE_CELERITY, SIMULATION_GRID_RESOLUTION, BOUNDARY_CONDITION, // Simulation parameters.
            {source}, // Source parameters.
            {obstacle}); // Obstacles parameters.

    app.Run(TIME_SCALE);
    return 0;
}