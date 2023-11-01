#ifndef MATERIAL
#define MATERIAL

#include "Vector.h"
#include "Random.h"
#include "JSON.h"
#include "Ray.h"
#include <iostream>

class Materials;

struct Record {
    double distance;
    Vector point, normal;
    Materials* objectMaterial;
    bool frontFace;

    void setNormal(const Ray& r, const Vector& normal) {
        if (r.direction.dot(normal) < 0) {
            frontFace = true;
            this -> normal = normal;
        } else {
            frontFace = false;
            this -> normal = -normal;
        }
    } 
};

class Materials {
    public:
    Vector reflect(const Vector& dir, const Vector& normal) const {
        return (dir) - ((normal * dir.dot(normal)) * 2);
    }

    virtual Vector calcEmission() const = 0;
    virtual bool scatter(Record& record, const Ray& r, Vector& attenuation, Ray& scattered) const = 0;
    virtual Vector getColor() const = 0;
    virtual void createMaterial(const Types& json) = 0; 
};

class Material : public Materials {
    public:
    Vector color;
    double specularProbability, smoothness;

    Material(const Vector& c = Color(255, 255, 255), const double& smooth = 0, const double& specProb = 0) {
        this -> color = c;
        this -> smoothness = smooth;
        this -> specularProbability = specProb;
    }

    Vector calcEmission() const override {
        return BLACK;
    }

    bool scatter(Record& record, const Ray& r, Vector& attenuation, Ray& scattered) const override {
        attenuation = color;
        if (specularProbability < randDbl()) {
            scattered = Ray(record.point, record.normal + randomInSphere().unitVector());
            if (scattered.direction.nearZero())
                return scatter(record, r, attenuation, scattered);
        } else {
            scattered = Ray(record.point, interpolate(record.normal + randomInSphere().unitVector(), 
                                                      reflect(r.direction, record.normal), smoothness));
        }

        return true;
    }

    Vector getColor() const override {
        return color;
    }

    void createMaterial(const Types& json) {
        Types* temp = json.jsonValue.get("color");
        if (temp != nullptr)
            color = Color(*(temp -> vecValue));
        
        temp = json.jsonValue.get("smoothness");
        if (temp != nullptr) 
            smoothness = temp -> dblValue;
        
        temp = json.jsonValue.get("specular_probability");
        if (temp != nullptr)
            specularProbability = temp -> dblValue;
    }

};

class Light : public Materials {
    public:
    Vector emissionColor;
    double lightStrength;

    Light(const Vector& c = Color(255, 255, 255), const double& lightStrength = 4) {
        this -> emissionColor = c * lightStrength;
        this -> lightStrength = lightStrength;
    }

    Vector calcEmission() const override {
        return emissionColor;
    }

    bool scatter(Record& record, const Ray& r, Vector& attenuation, Ray& scattered) const override {
        return false;
    }

    Vector getColor() const override {
        return BLACK;
    }

    void createMaterial(const Types& json) {
        Light* material = new Light();
        Types* temp = json.jsonValue.get("brightness");
        if (temp != nullptr)
            lightStrength = temp -> dblValue;
        
        temp = json.jsonValue.get("emission_color");
        if (temp != nullptr) 
            emissionColor = Color(*(temp -> vecValue)) * lightStrength;
    }
};


Materials* buildMaterial(const Types& json) {
    Types* matType = json.jsonValue.get("material_type");
    
    if (matType == nullptr || matType -> strValue == "material") { //Solid or Metal - Default
        Material* material = new Material();
        material -> createMaterial(json.jsonValue);
        return material;
    } else { //Light
        Light* material = new Light();
        material -> createMaterial(json.jsonValue);
        return material;
    }
}

#endif