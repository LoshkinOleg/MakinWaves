#pragma once

#include <functional>
#include <cassert>
#include <iostream>
#include <algorithm>

struct Vec2f {
    Vec2f() = default;
    Vec2f(const double x, const double y) : x(x), y(y) {}

    double x = 0, y = 0;
};

class Rect {
public:
    Rect() = delete;
    Rect(const Vec2f min, const Vec2f max) : min_(min), max_(max) {}

    inline Vec2f GetMin() const {
        return min_;
    }
    inline Vec2f GetMax() const {
        return max_;
    }

private:
    Vec2f min_ = {}, max_ = {};
};

typedef double SeekPosition;
typedef double Sample;
typedef std::function<Sample(SeekPosition)> SourceCallback;

class Source {
public:
    Source() = delete;
    Source(const SourceCallback callback, const Vec2f pos) : callback(callback), pos_(pos) {}

    void SetPos(const Vec2f newPos) {
        pos_ = newPos;
    }

    inline Vec2f GetPos() const {
        return pos_;
    }

    inline Sample GetSample(const SeekPosition seek) const {
        return callback(seek);
    }

    const SourceCallback callback;
private:
    Vec2f pos_ = {};
};

typedef std::vector<std::vector<double>> DisplacementField2D;

enum BoundaryCondition {
    DIRICHLET = 0, // Displacement at boundaries is always 0. Causes perfect reflections (no loss of energy).
    MUR = 1 // Boundaries absorb energy. Absorbs most of the energy (but not all).
};
enum FieldUpdate {
    HAROON = 0, // Taken from Haroon Stephen's video course: https://youtu.be/O6fqBxuM-g8
};

// Implemented following Haroon Stephen's videos: https://youtu.be/MqlDQzw3x6Q
class ExplicitWorld2D {
public:
    ExplicitWorld2D() = delete;
    ExplicitWorld2D(
        const double celerity,
        const DisplacementField2D& initialConditions,
        const FieldUpdate fieldUpdateType,
        const BoundaryCondition boundaryConditionType,
        const std::vector<Source>& sources = {},
        const std::vector<Rect>& rects = {}
    ) : resolutionX(initialConditions.size()), resolutionY(initialConditions[0].size()),
        celerity(celerity),
        dLast_(initialConditions), dCurrent_(initialConditions), dNext_(initialConditions),
        rects_(rects), sources_(sources),
        boundaryConditionType(boundaryConditionType),
        fieldUpdateType(fieldUpdateType)
    {
        assert(initialConditions.size() > 0 && initialConditions[0].size() > 0 && "The initial displacement field should have at least a 2x2 size.");
    }

    /*
    @return: current displacement field.
    */
    const DisplacementField2D& Update(const double dt)
    {
        // Courant–Friedrichs–Lewy condition: https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
        const Vec2f cfl =
        {
            (celerity * dt) / (1.0f / (double)resolutionX) ,
            (celerity * dt) / (1.0f / (double)resolutionY)
        };
        if (cfl.x + cfl.y > 1.0f)
        {
            // std::cout << "Courant-Freidrichs-Lewy condition is not fulfilled, simulation is inaccurate, skipping...\n";
            // return dNext_;
        }

        /*
            Neighboring cells:
                   |  X - 1  |   X     | X + 1
             ----- | --------------------------
             Y - 1 | -       | Top     | -
                 Y | Left    | Center  | Right
             Y + 1 | -       | Bottom  | -
        */

        // Update the bulk of the field.
        switch (fieldUpdateType)
        {
        case FieldUpdate::HAROON:
        {
            for (size_t x = 1; x < resolutionX - 1; x++)
            {
                for (size_t y = 1; y < resolutionY - 1; y++)
                {
                    const double currentCenter = dCurrent_[x][y];
                    const double lastCenter = dLast_[x][y];
                    const double currentRight = dCurrent_[x + 1][y];
                    const double currentLeft = dCurrent_[x - 1][y];
                    const double currentTop = dCurrent_[x][y + 1];
                    const double currentBottom = dCurrent_[x][y - 1];

                    // Formula: dNextCenter = 2*dCurrentCenter - dLastCenter + CFL^2 * (dCurrentRight + dCurrentLeft - 4*dCurrentCenter + dCurrentBottom + dCurrentTop) + (dt^2 * source@t)
                    dNext_[x][y] =
                        2.0f * currentCenter
                        - lastCenter
                        + ((cfl.x + cfl.y) * (cfl.x + cfl.y))
                        * (currentRight + currentLeft - 4.0f * currentCenter + currentTop + currentBottom);
                }
            }
        }break;
        default:break;
        }

        // Update field with forcing functions (aka sound sources). Using Haroon's discretization: https://www.youtube.com/watch?v=O6fqBxuM-g8
        for (auto& src : sources_)
        {
            const size_t pixelPosX = std::clamp(src.GetPos().x, 0.0, 1.0) * resolutionX;
            const size_t pixelPosY = std::clamp(src.GetPos().y, 0.0, 1.0) * resolutionY;

            // Formula: dNextCenter += dt^2 * f(t).
            // dNext_[pixelPosX][pixelPosY] += dt * dt * src.GetSample(t_);
            dNext_[pixelPosX][pixelPosY] = src.GetSample(t_);
        }

        // Update boundary conditions.
        switch (boundaryConditionType)
        {
        case BoundaryCondition::DIRICHLET:
        {
            // Note: unless initial conditions at boundaries are not at 0,
            // this bit of code should never actually chage anything since
            // the boundaries are not updated in the code above.
            // And if initial conditions at boundaries are not 0, this code
            // should only actually change the displacement values only once.

            // Left boundary. x = left boundary, y varies.
            for (size_t y = 0; y < resolutionY; y++)
            {
                // Formula: dNextCenter = 0
                dNext_[0][y] = 0.0f;
            }
            // Right boundary. x = right boundary, y varies.
            for (size_t y = 0; y < resolutionY; y++)
            {
                // Formula: dNextCenter = 0
                dNext_[resolutionX - 1][y] = 0.0f;
            }
            // Top boundary. x varies, y = top boundary.
            for (size_t x = 0; x < resolutionX; x++)
            {
                // Formula: dNextCenter = 0
                dNext_[x][0] = 0.0f;
            }
            // Bottom boundary. x varies, y = bottom boundary.
            for (size_t x = 0; x < resolutionX; x++)
            {
                // Formula: dNextCenter = 0
                dNext_[x][resolutionY - 1] = 0.0f;
            }

            // Obstacles.
            for (auto& obj : rects_)
            {
                const size_t minPixelPosX = std::clamp(obj.GetMin().x, 0.0, 1.0) * resolutionX;
                const size_t minPixelPosY = std::clamp(obj.GetMin().y, 0.0, 1.0) * resolutionY;
                const size_t maxPixelPosX = std::clamp(obj.GetMax().x, 0.0, 1.0) * resolutionX;
                const size_t maxPixelPosY = std::clamp(obj.GetMax().y, 0.0, 1.0) * resolutionY;

                for (size_t x = minPixelPosX; x < maxPixelPosX; x++)
                {
                    for (size_t y = minPixelPosY; y < maxPixelPosY; y++)
                    {
                        // Formula: dNextCenter = 0
                        dNext_[x][y] = 0.0f;
                    }
                }
            }
        } break;
        case BoundaryCondition::MUR:
        {
            // Taken from: https://youtu.be/MqlDQzw3x6Q

            // Left boundary. x = left boundary, y varies.
            for (size_t y = 0; y < resolutionY; y++)
            {
                const double currentRight = dCurrent_[1][y];
                const double center = dCurrent_[0][y];
                const double nextRight = dNext_[1][y];

                // Formula: dNextCenter = dCurrentRight + (CFL - 1 / CFL + 1) * (dNextRight - dCurrentCenter)
                dNext_[0][y] =
                    currentRight
                    + (((cfl.x + cfl.y) - 1.0f) / ((cfl.x + cfl.y) + 1.0f))
                    * (nextRight - center);
            }
            // Right boundary. x = right boundary, y varies.
            for (size_t y = 0; y < resolutionY; y++)
            {
                const double currentLeft = dCurrent_[resolutionX - 2][y];
                const double center = dCurrent_[resolutionX - 1][y];
                const double nextLeft = dNext_[resolutionX - 2][y];

                // Formula: dNextCenter = dCurrentLeft + (CFL - 1 / CFL + 1) * (dNextLeft - dCurrentCenter)
                dNext_[resolutionX - 1][y] =
                    currentLeft
                    + (((cfl.x + cfl.y) - 1.0f) / ((cfl.x + cfl.y) + 1.0f))
                    * (nextLeft - center);
            }
            // Top boundary. x varies, y = top boundary.
            for (size_t x = 0; x < resolutionX; x++)
            {
                const double currentBottom = dCurrent_[x][1];
                const double current = dCurrent_[x][0];
                const double nextBottom = dNext_[x][1];

                // Formula: dNextCenter = dCurrentBottom + (CFL - 1 / CFL + 1) * (dNextBottom - dCurrentCenter)
                dNext_[x][0] =
                    currentBottom
                    + (((cfl.x + cfl.y) - 1.0f) / ((cfl.x + cfl.y) + 1.0f))
                    * (nextBottom - current);
            }
            // Bottom boundary. x varies, y = bottom boundary.
            for (size_t x = 0; x < resolutionX; x++)
            {
                const double currentTop = dCurrent_[x][resolutionY - 2];
                const double current = dCurrent_[x][resolutionY - 1];
                const double nextTop = dNext_[x][resolutionY - 2];

                // Formula: dNextCenter = dCurrentTop + (CFL - 1 / CFL + 1) * (dNextTop - dCurrentCenter)
                dNext_[x][resolutionY - 1] =
                    currentTop
                    + (((cfl.x + cfl.y) - 1.0f) / ((cfl.x + cfl.y) + 1.0f))
                    * (nextTop - current);
            }

            // Obstacles.
            for (auto& obj : rects_)
            {
                const size_t minPixelPosX = (size_t)std::ceilf(std::clamp(obj.GetMin().x, 0.0, 1.0) * (double)resolutionX);
                const size_t minPixelPosY = (size_t)std::ceilf(std::clamp(obj.GetMin().y, 0.0, 1.0) * (double)resolutionY);
                const size_t maxPixelPosX = (size_t)std::ceilf(std::clamp(obj.GetMax().x, 0.0, 1.0) * (double)resolutionX);
                const size_t maxPixelPosY = (size_t)std::ceilf(std::clamp(obj.GetMax().y, 0.0, 1.0) * (double)resolutionY);

                // Using inverted Mur for the boundaries of the obstacle.

                // Left boundary. x = left boundary, y varies.
                for (size_t y = minPixelPosY; y < maxPixelPosY; y++)
                {
                    const double currentLeft = dCurrent_[minPixelPosX - 1][y];
                    const double center = dCurrent_[minPixelPosX][y];
                    const double nextLeft = dNext_[minPixelPosX - 1][y];

                    // Formula: dNextCenter = dCurrentLeft + (CFL - 1 / CFL + 1) * (dNextLeft - dCurrentCenter)
                    dNext_[minPixelPosX][y] =
                        currentLeft
                        + (((cfl.x + cfl.y) - 1.0f) / ((cfl.x + cfl.y) + 1.0f))
                        * (nextLeft - center);
                }
                // Right boundary. x = right boundary, y varies.
                for (size_t y = minPixelPosY; y < maxPixelPosY; y++)
                {
                    const double currentRight = dCurrent_[maxPixelPosX + 1][y];
                    const double center = dCurrent_[maxPixelPosX][y];
                    const double nextRight = dNext_[maxPixelPosX + 1][y];

                    // Formula: dNextCenter = dCurrentRight + (CFL - 1 / CFL + 1) * (dNextRight - dCurrentCenter)
                    dNext_[maxPixelPosX][y] =
                        currentRight
                        + (((cfl.x + cfl.y) - 1.0f) / ((cfl.x + cfl.y) + 1.0f))
                        * (nextRight - center);
                }
                // Top boundary. x varies, y = top boundary.
                for (size_t x = minPixelPosX; x < maxPixelPosX; x++)
                {
                    const double currentTop = dCurrent_[x][minPixelPosY - 1];
                    const double center = dCurrent_[x][minPixelPosY];
                    const double nextTop = dNext_[x][minPixelPosY - 1];

                    // Formula: dNextCenter = dCurrentTop + (CFL - 1 / CFL + 1) * (dNextTop - dCurrentCenter)
                    dNext_[x][minPixelPosY] =
                        currentTop
                        + (((cfl.x + cfl.y) - 1.0f) / ((cfl.x + cfl.y) + 1.0f))
                        * (nextTop - center);
                }
                // Bottom boundary. x varies, y = bottom boundary.
                for (size_t x = minPixelPosX; x < maxPixelPosX; x++)
                {
                    const double currentBottom = dCurrent_[x][maxPixelPosY + 1];
                    const double center = dCurrent_[x][maxPixelPosY];
                    const double nextBottom = dNext_[x][maxPixelPosY + 1];

                    // Formula: dNextCenter = dCurrentBottom + (CFL - 1 / CFL + 1) * (dNextBottom - dCurrentCenter)
                    dNext_[x][maxPixelPosY] =
                        currentBottom
                        + (((cfl.x + cfl.y) - 1.0f) / ((cfl.x + cfl.y) + 1.0f))
                        * (nextBottom - center);
                }

                // Using Dirichlet for the unreachable area of the obstacle.
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

        // Advance displacement fields buffers in time.
        t_ += dt;
        std::copy(dCurrent_.begin(), dCurrent_.end(), dLast_.begin());
        std::copy(dNext_.begin(), dNext_.end(), dCurrent_.begin());
        return dNext_;
    }

    const size_t resolutionX, resolutionY;
    const double celerity;
    const BoundaryCondition boundaryConditionType;
    const FieldUpdate fieldUpdateType;
private:
    double t_ = 0.0f;
    std::vector<Rect> rects_ = {};
    std::vector<Source> sources_ = {};
    DisplacementField2D dLast_ = {}, dCurrent_ = {}, dNext_ = {}; // 3 frames of displacement.
};