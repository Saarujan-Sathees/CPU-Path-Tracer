#ifndef RAY
#define RAY

#include "Vector.h"

struct Ray {
    Vector origin, direction;

    Ray() {
        origin = Vector();
        direction = Vector();
    }

    Ray(const Vector& origin, const Vector& direction) {
        this -> origin = origin;
        this -> direction = direction;
    }

    inline Vector point(const double& t) const {
        return origin + direction * t;
    }
};

#endif