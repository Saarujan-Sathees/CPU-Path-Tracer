#ifndef RENDER
#define RENDER

#include "Vector.h"
#include "Objects.h"
#include "Materials.h"
#include "Scene.h"
#include <iostream>

Vector traceRay(const Ray& r, const Scene& scene, const Vector& backgroundColor = BLACK, const short& bounceCount = 50) {
    Record record;
    if (bounceCount <= 0) return BLACK;
    if (!scene.intersect(r, 0.001, INFINITY, record)) return backgroundColor;

    Vector emitted = record.objectMaterial -> calcEmission(), color;
    Ray scattered;
    
    if (!record.objectMaterial -> scatter(record, r, color, scattered))
        return emitted;
    
    return emitted + color * traceRay(scattered, scene, backgroundColor, bounceCount - 1);
    /*if (scene.intersect(r, 0.001, INFINITY, record))
        return traceRay(Ray(record.point, record.normal + randomInSphere().unitVector()), scene, bounceCount - 1) * 0.5;
    else {
        Vector unit_direction = r.direction.unitVector();
        auto t = 0.5 * (unit_direction[1] + 1.0);
        return Vector(1.0 - t, 1.0 - t, 1.0 - t) + Vector(0.5, 0.7, 1.0) * t;
    }*/
}

#endif