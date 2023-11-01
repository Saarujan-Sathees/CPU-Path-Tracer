#ifndef RANDOM
#define RANDOM

#include <random>
#include "Vector.h"

static std::minstd_rand randEngine(std::random_device{}());
static std::uniform_real_distribution<double> zeroToOneDist(0, 1);
static std::uniform_real_distribution<double> oneToOneDist(-1, 1);

inline double randDbl(double min = 0, double max = 1.0) {
    if (min == 0 && max == 1)
        return zeroToOneDist(randEngine);
    else if (min == -1 && max == 1)
        return oneToOneDist(randEngine);

    static std::uniform_real_distribution<double> distribution(min, max);
    return distribution(randEngine);
}

inline Vector randVector(double min = 0, double max = 1) {
    return Vector(randDbl(min, max), randDbl(min, max), randDbl(min, max));
}

inline Vector randomInSphere() {
    Vector point;
    do {
        point = randVector(-1, 1);
    } while (point.squaredLength() >= 1);

    return point;
}

inline Vector randomInCircle() {
    Vector point;
    do {
        point = Vector(randDbl(-1, 1), randDbl(-1, 1), 0);
    } while (point.squaredLength() >= 1);
    
    return point;
}

#endif