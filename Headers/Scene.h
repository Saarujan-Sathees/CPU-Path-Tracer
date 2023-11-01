#ifndef SCENE
#define SCENE

#include "Objects.h"
#include "Materials.h"
#include <vector>

class Scene : public Traceable {
    private:
    std::vector<Traceable*> objects;

    public:
    Scene() {
        objects = std::vector<Traceable*>(0);
    }

    Scene(Traceable* o) {
        objects.emplace_back(o);
    }

    void push(Traceable* o) {
        objects.emplace_back(o);
    }

    Traceable* pop() {
        Traceable* res = objects[objects.size() - 1];
        objects.pop_back();
        return res;
    }

    bool intersect(const Ray& r, const double& min, const double& max, Record& info) const {
        Record temp;
        double closestDistance = max;
        bool hitObject = false;
        for (int i = 0; i < size(); ++i) {
            if (objects[i] -> intersect(r, min, closestDistance, temp)) {
                closestDistance = temp.distance;
                hitObject = true;
                info = temp;
            }
        }

        return hitObject;
    }

    int size() const {
        return objects.size();
    }

    void buildObject(const Types& json) {
        return;
    }

    void printType() const {
        std::cout << "Scene\n";
    }

};

#endif