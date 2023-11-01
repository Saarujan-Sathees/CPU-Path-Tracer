#ifndef CAMERA
#define CAMERA

#include "Vector.h"
#include "Ray.h"
#include "Random.h"
#include <math.h>

struct Camera {
    private:
    static const double degToHalfRad;

    public: 
    Vector origin, target, horizontal, vertical, adjustedLLC, u, v;
    short fieldOfView;
    double lensRadius, focusDist;

    Camera() {}

    Camera(const Vector& o, const Vector& tP, const short& FOV, const double& aspectRatio, const double& fD = -1, const double& a = 0) {        
        origin = o;
        target = tP;
        fieldOfView = FOV;
        lensRadius = a * 0.5;
        if (fD == -1) 
            focusDist = (origin - tP).length();
        else 
            focusDist = fD;

        Vector w = (origin - tP).unitVector();
        u = Vector(0, 1, 0).cross(w);
        v = u.cross(w);

        double viewportHeight = 2 * std::tan(FOV * degToHalfRad);
        horizontal = u * focusDist * viewportHeight * aspectRatio;
        vertical = v * focusDist * viewportHeight;
        adjustedLLC = horizontal * -0.5 - vertical * 0.5 - w * focusDist;
    }

    inline Ray generateRay(const double& x, const double& y) const {
        Vector rand = randomInCircle() * lensRadius;
        Vector offset = u * rand[0] + v * rand[1];

        return Ray(origin + offset, adjustedLLC + horizontal * x + vertical * y - offset);
    }
};

const double Camera::degToHalfRad = 3.1415926535 / 360;

#endif