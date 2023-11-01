#ifndef OBJECTS
#define OBJECTS

#include "Vector.h"
#include "Ray.h"
#include "Materials.h"
#include "JSON.h"

#define EPSILON 1E-16

class Traceable {
    public:
    virtual bool intersect(const Ray& r, const double& min, const double& max, Record& info) const = 0;
    virtual void buildObject(const Types& json) = 0;
    virtual void printType() const = 0;
};

class BoundingBox {
    public:
    Vector minimum, maximum;

    BoundingBox() {}
    BoundingBox(const Vector& minimum, const Vector& maximum) {
        this -> minimum = minimum;
        this -> maximum = maximum;
    }

    bool intersect(const Ray& r, const double& min, const double& max) const {
        double recDirection, tMin = min, tMax = max;
        for (int i = 0; i < 3; ++i) {
            recDirection = 1.0 / r.direction[i];
            if (recDirection < 0.0) {
                tMin = std::max((maximum[i] - r.origin[i]) * recDirection, tMin);
                tMax = std::min((minimum[i] - r.origin[i]) * recDirection, tMax);
            } else {
                tMin = std::max((minimum[i] - r.origin[i]) * recDirection, tMin);
                tMax = std::min((maximum[i] - r.origin[i]) * recDirection, tMax);
            }
            
            if (tMax <= tMin)
                return false;
        }

        return true;
    }

    void printType() const {
        std::cout << "Bounding Box\n";
    }

};

class Sphere : public Traceable {
    public:
    Materials* material;
    Vector center;
    BoundingBox* boundary;
    double radius, radSqr, recRad;

    Sphere() {
        material = new Material();
        center = Vector();
        radius = 0;
        radSqr = 0;
        recRad = 0;
    }

    Sphere(const Vector& center, const double& radius, Materials* material) {
        this -> center = center;
        this -> radius = radius;
        this -> material = material;
        radSqr = radius * radius;
        recRad = 1.0 / radius;
        boundary = new BoundingBox(center - radius, center + radius);
    }

    inline bool intersect(const Ray& r, const double& min, const double& max, Record& info) const {
        if (!boundary -> intersect(r, min, max)) return false;

        const Vector adjustedCenter = r.origin - center;
        double a = r.direction.squaredLength();
        double b = r.direction.dot(adjustedCenter);

        double disc = b * b - a * (adjustedCenter.squaredLength() - radSqr);
        if (disc < 0) //No Solutions
            return false;
        
        a = 1.0 / a; //Storing the reciprocal to reduce the amount of divisions
        disc = std::sqrt(disc); //Storing the root of discriminant to reduce the amount of root calculations
        double distance = (-b - disc) * a;
        if (distance < min || distance > max) {
            distance = (-b + disc) * a;
            if (distance < min || distance > max) return false;
        }
            
        info.distance = distance;
        info.point = r.point(distance);
        info.objectMaterial = material;
        info.setNormal(r, (info.point - center) * recRad);
        return true;
    }

    void buildObject(const Types& json) {
        Types* temp = json.jsonValue.get("position");
        if (temp != nullptr) 
            center = *(temp -> vecValue);
        else
            center = Vector(0, 0, 0);
        
        temp = json.jsonValue.get("radius");
        if (temp != nullptr) 
            radius = temp -> dblValue;
        else
            radius = 1;
        
        radSqr = radius * radius;
        recRad = 1.0 / radius;
        boundary = new BoundingBox(center - radius, center + radius);
        material = buildMaterial(json);
    }

    void printType() const {
        std::cout << "Sphere\n";
    }

};

class Triangle : public Traceable {
    public:
    Materials* material;
    Vector vOne, vTwo, vThree, planeNormal, a, b;

    Triangle() {}
    Triangle(const Vector& v1, const Vector& v2, const Vector& v3, Materials* material) {
        vOne = v1;
        vTwo = v2;
        vThree = v3;
        this -> material = material;

        a = vTwo - vOne;
        b = vThree - vOne;
        planeNormal = a.cross(b);
    }

    inline bool intersect(const Ray& r, const double& min, const double& max, Record& info) const {
        double recDet = -r.direction.dot(planeNormal);
        if (recDet < EPSILON) 
            return false;

        recDet = 1.0 / recDet;
        Vector adjustedOrigin  = r.origin - vOne, detAO = adjustedOrigin.cross(r.direction);
        double u = b.dot(detAO) * recDet;
        double v = -a.dot(detAO) * recDet, w = 1 - u - v;
        info.distance = adjustedOrigin.dot(planeNormal) * recDet;

        if (u < 0 || v < 0 || w < 0 || info.distance < min || info.distance > max) 
            return false;

        info.objectMaterial = material;
        info.point = r.point(info.distance);
        info.setNormal(r, planeNormal.unitVector());
        
        return true;
    }

    void correctNormal(const Vector& cameraOrigin) {
        Vector rayDirection = vOne - cameraOrigin;

        if (-rayDirection.dot(planeNormal) < EPSILON) {
            Vector temp = vTwo;
            vTwo = vThree;
            vThree = temp;

            a = vTwo - vOne;
            b = vThree - vOne;
            planeNormal *= -1;
        }
    }

    void buildObject(const Types& json) {
        Types* temp = json.jsonValue.get("vertices");
        if (temp != nullptr) {
            vOne = *(temp -> arrValue[0] -> vecValue);
            vTwo = *(temp -> arrValue[1] -> vecValue);
            vThree = *(temp -> arrValue[2] -> vecValue);
        } else {
            vOne = Vector(3, 6, 0);
            vTwo = Vector(6, 0, 0);
            vThree = Vector(0, 0, 0);
        }
        
        a = vTwo - vOne;
        b = vThree - vOne;
        planeNormal = a.cross(b);
        material = buildMaterial(json);
    }

    void printType() const {
        std::cout << "Triangle\n";
    }

};

class Quadrilateral : public Traceable {
    public:
    Triangle* triOne, *triTwo;

    Quadrilateral() {}
    Quadrilateral(const Vector& v1, const Vector& v2, const Vector& v3, const Vector& v4, Materials* material) {
        double distTwo = (v1 - v2).squaredLength(), distThree = (v1 - v3).squaredLength(), distFour = (v1 - v4).squaredLength();
        double largestDist = std::max(distTwo, std::max(distThree, distFour));

        if (largestDist == distTwo) {
            triOne = new Triangle(v1, v3, v4, material);
            triTwo = new Triangle(v4, v3, v2, material);
        } else if (largestDist == distThree) {
            triOne = new Triangle(v1, v2, v4, material);
            triTwo = new Triangle(v4, v3, v2, material);
        } else {
            triOne = new Triangle(v1, v2, v3, material);
            triTwo = new Triangle(v2, v3, v4, material);
        }
    }

    inline bool intersect(const Ray& r, const double& min, const double& max, Record& info) const {
        return triOne -> intersect(r, min, max, info) || triTwo -> intersect(r, min, max, info);
    }

    void correctNormal(const Vector& cameraOrigin) {
        triOne -> correctNormal(cameraOrigin);
        triTwo -> correctNormal(cameraOrigin);
    }

    void buildObject(const Types& json) {
        Types* temp = json.jsonValue.get("vertices");
        Materials* material = buildMaterial(json);

        if (temp != nullptr) {
            Vector v1 = *(temp -> arrValue[0] -> vecValue), v2 = *(temp -> arrValue[1] -> vecValue),
                   v3 = *(temp -> arrValue[2] -> vecValue), v4 = *(temp -> arrValue[3] -> vecValue);
            
            double distTwo = (v1 - v2).squaredLength(), distThree = (v1 - v3).squaredLength(), distFour = (v1 - v4).squaredLength();
            double largestDist = std::max(distTwo, std::max(distThree, distFour));

            if (largestDist == distTwo) {
                triOne = new Triangle(v1, v3, v4, material);
                triTwo = new Triangle(v4, v3, v2, material);
            } else if (largestDist == distThree) {
                triOne = new Triangle(v1, v2, v4, material);
                triTwo = new Triangle(v4, v3, v2, material);
            } else {
                triOne = new Triangle(v1, v2, v3, material);
                triTwo = new Triangle(v2, v3, v4, material);
            }
        } else {
            triOne = new Triangle(Vector(0, 0, 0), Vector(0, 6, 0), Vector(6, 6, 0), material);
            triTwo = new Triangle(Vector(0, 0, 0), Vector(0, 6, 0), Vector(6, 0, 0), material);
        }
    }

    void printType() const {
        std::cout << "Quadrilateral\n";
    }

};

class AlignedRectangle : public Traceable {
    private:
    double aOne, aTwo, bOne, bTwo, k, axis, aTO, bTO, area;
    int aAx = 0, bAx = 1;
    Vector normal = Vector(0, 0, 1);
    Materials* matType;

    public:
    AlignedRectangle() {}
    AlignedRectangle(double xO, double xT, double yO, double yT, double z, int ax, Materials* mat) {
        aOne = xO;
        aTwo = xT;
        bOne = yO;
        bTwo = yT;
        k = z;
        axis = ax;
        aTO = aTwo - aOne;
        bTO = bTwo - bOne;
        area = aTO * bTO;
        if (axis == 0) {
            aAx = 1;
            bAx = 2;
            normal = Vector(1, 0, 0);
        } else if (axis == 1) {
            bAx = 2;
            normal = Vector(0, 1, 0);
        }
        matType = mat;
    }

    bool intersect(const Ray& r, const double& min, const double& max, Record& info) const {
        const double t = (k - r.origin[axis]) / r.direction[axis];
        if (t < min || t > max) return false;

        const double a = r.origin[aAx] + t * r.direction[aAx];
        const double b = r.origin[bAx] + t * r.direction[bAx];
        if (a < aOne || a > aTwo || b < bOne || b > bTwo) return false;
        
        info.distance = t;
        info.point = r.point(info.distance);
        info.setNormal(r, normal);
        info.objectMaterial = matType;
        return true;
    }

    void buildObject(const Types& json) {
        Types* temp = json.jsonValue.get("values");
        if (temp != nullptr) {
            aOne = temp -> arrValue[0] -> dblValue;
            aTwo = temp -> arrValue[1] -> dblValue;

            bOne = temp -> arrValue[2] -> dblValue;
            bTwo = temp -> arrValue[3] -> dblValue;
        } 
        
        temp = json.jsonValue.get("k");
        if (temp != nullptr)
            k = temp -> dblValue;
        
        temp = json.jsonValue.get("axis");
        if (temp != nullptr)
            axis = temp -> dblValue;
        
        matType = buildMaterial(json);
        aTO = aTwo - aOne;
        bTO = bTwo - bOne;
        area = aTO * bTO;

        if (axis == 0) {
            aAx = 1;
            bAx = 2;
            normal = Vector(1, 0, 0);
        } else if (axis == 1) {
            bAx = 2;
            normal = Vector(0, 1, 0);
        }
    }

    void printType() const {
        std::cout << "Axis-Aligned Rectangle\n";
    }

};

#endif