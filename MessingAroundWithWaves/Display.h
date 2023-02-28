#pragma once

#include <cmath>
#include <vector>
#include <functional>

#include "raylib.h"

inline double ContrastExponent(const double x, const double exponent)
{
    return 1.0f / (1.0f + std::exp(-exponent * x)) - 0.5f;
}

inline double ContrastReinhardt(const double x, const double t)
{
    return x * t / (x * (t - 1.0) + 1.0);
}

Color ToColor(const double displacement) {

    if (displacement > 1.0f) {
        return
        {
            255,
            255,
            255,
            255
        };
    }
    else if (displacement < 0.0f) {
        return { 0,0,0,255 };
    }
    else
    {
        return
        {
            (unsigned char)(displacement * 255.0f),
            (unsigned char)(displacement * 255.0f),
            (unsigned char)(displacement * 255.0f),
            255
        };
    }

}

void DrawDisplacement(const std::vector<std::vector<double>>& d, const size_t displayWidth, const size_t displayHeight) {
    for (size_t x = 0; x < displayWidth; x++)
    {
        for (size_t y = 0; y < displayHeight; y++)
        {
            const double fracX = x / (float)displayWidth;
            const double fracY = y / (float)displayHeight;
            const size_t coordX = d.size() * fracX;
            const size_t coordY = d[0].size() * fracY;
            DrawPixel(x, y, ToColor(ContrastReinhardt(d[coordX][coordY], 10)));
        }
    }
    
    // for (size_t x = 0; x < d.size(); x++)
    // {
    //     for (size_t y = 0; y < d[0].size(); y++)
    //     {
    //         DrawPixel(x, y, ToColor(d[x][y], contrastExponent));
    //     }
    // }
}