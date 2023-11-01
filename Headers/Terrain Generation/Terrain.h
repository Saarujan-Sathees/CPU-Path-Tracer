#ifndef TERRAIN
#define TERRAIN

#include "../Vector.h"
#include "../Objects.h"
#include "../JSON.h"
#include "LookupTable.h"
#include "Noise.h"

const Vector xOff = Vector(1, 0, 0), yOff = Vector(0, 1, 0), zOff = Vector(0, 0, 1);

class Terrain : public Traceable {
    Vector Dimensions, cornerOffset;
    std::vector<std::vector<Triangle*>> Faces;
    std::vector<BoundingBox*> boundaries;
    BoundingBox* mainBoundary;
    SimplexNoise* noiseMap;
    Materials* material;
    std::vector<Vector> Coord;
    short cubesPerChunk;
    double maxHeight, halfHeight;

    inline double Density(const Vector &c) {
        Vector bound = c;
        if (c[0] >= Dimensions[0] || c[1] >= Dimensions[1] || c[2] >= Dimensions[2]) 
            bound = Dimensions - Vector(1, 1, 1);
        else if (c[0] < 0 || c[1] < 0 || c[2] < 0) 
            bound = Vector(0, 0, 0);
  
        double density = -bound[1] + halfHeight + noiseMap -> valueAt(bound);

        return density;
    }

    inline Vector calcNormal(const Vector &c) {
        Vector bound = c;
        if (c[0] >= Dimensions[0] - 1 || c[1] >= Dimensions[1] || c[2] >= Dimensions[2] - 1) 
            bound = Dimensions - Vector(1, 1, 1);
        if (c[0] < 0 || c[1] < 0 || c[2] < 0) 
            bound = Vector(0, 0, 0);

        Vector norm(Density(bound - xOff) - Density(bound + xOff), 
                   Density(bound - yOff) - Density(bound + yOff),
                   Density(bound - zOff) - Density(bound + zOff));

        return norm / std::sqrt(norm.dot(norm));
    }

    inline Vector Lerp(const Vector &cA, const Vector &cB, const double &t) const {
        return cA + (cB - cA) * t;
    }

    inline Vector Interpolate(const int &cA, const int &cB) {
        const Vector posA = Coord[cA], posB = Coord[cB];
        const double densityA = Density(posA), t = (maxHeight - densityA) / (Density(posB) - densityA);
        //const Vector normA = calcNormal(posA), normB = calcNormal(posB);
        
        //Vertex v;
        //v.vXYZ = Lerp(toWorld(posA), toWorld(posB), t);
        //v.tNormal = Lerp(normA, normB, t);
        //v.tNormal = v.tNormal / std::sqrt(v.tNormal.dot(v.tNormal));
        return Lerp(posA, posB, t); 
    }

    void GenerateChunk(const Vector &chunkCoord, const Vector &cubesInChunk, const int& chunkIndex, const Vector& camOrigin) {
        char edgeOne, edgeTwo, edgeThr;
        double minY = INFINITY, maxY = -1E16;
        std::vector<char> l;
        Triangle* temp = new Triangle();
        Vector cOne, cTwo, cThree, cFour;
                
        for (short x = chunkCoord[0]; x < cubesInChunk[0] + chunkCoord[0]; ++x) {
            for (short y = chunkCoord[1]; y < cubesInChunk[1] + chunkCoord[1]; ++y) {
                /*cOne = Vector(x, 0, z), cTwo = Vector(x + 1, 0, z), cThree = Vector(x, 0, z + 1), cFour = Vector(x + 1, 0, z + 1);
                cOne.set(1, noiseMap -> valueAt(cOne[0], cOne[2]));
                cTwo.set(1, noiseMap -> valueAt(cTwo[0], cTwo[2]));
                cThree.set(1, noiseMap -> valueAt(cThree[0], cThree[2]));
                cFour.set(1, noiseMap -> valueAt(cFour[0], cFour[2]));

                temp = new Triangle(cOne + cornerOffset, cTwo + cornerOffset, cThree + cornerOffset, material);
                temp -> correctNormal(camOrigin);
                Faces[chunkIndex].push_back(temp);

                temp = new Triangle(cTwo + cornerOffset, cThree + cornerOffset, cFour + cornerOffset, material);
                temp -> correctNormal(camOrigin);
                Faces[chunkIndex].push_back(temp);

                minY = std::min(cOne[1], std::min(cTwo[1], std::min(cThree[1], std::min(cFour[1], minY))));
                maxY = std::max(cOne[1], std::max(cTwo[1], std::max(cThree[1], std::max(cFour[1], maxY))));*/
                
                for (short z = chunkCoord[2]; z < cubesInChunk[2] + chunkCoord[2]; ++z) {
                    Vector c = Vector(x, y, z) + chunkCoord;
                    Coord = { c, c + Vector(1, 0, 0), c + Vector(1, 0, 1), c + Vector(0, 0, 1), 
                              c + Vector(0, 1, 0), c + Vector(1, 1, 0), 
                              c + 1, c + Vector(0, 1, 1) };
                    
                    short Index = 0;
                    for (short i = 0; i < 8; ++i) {
                        if (Density(Coord[i]) < maxHeight) {
                            Index |= (1 << i);
                        }
                    }

                    if (Index == 0 || Index == 255) continue;

                    l = lookup[Index];
                    for (short i = 0; i < l.size(); i += 3) {
                        edgeOne = l[i], edgeTwo = l[i + 1], edgeThr = l[i + 2];
                        if (chunkCoord[0] == 0 && chunkCoord[2] == 0 && i == 0) //Calculating corner offset
                            cornerOffset -= Interpolate(cornerFromA[edgeOne], cornerFromB[edgeOne]);

                        temp = new Triangle(Interpolate(cornerFromA[edgeOne], cornerFromB[edgeOne]) + cornerOffset,
                                            Interpolate(cornerFromA[edgeTwo], cornerFromB[edgeTwo]) + cornerOffset,
                                            Interpolate(cornerFromA[edgeThr], cornerFromB[edgeThr]) + cornerOffset, material);

                        temp -> correctNormal(camOrigin);
                        minY = std::min(temp -> vOne[1], std::min(temp -> vTwo[1], std::min(temp -> vThree[1], minY)));
                        maxY = std::max(temp -> vOne[1], std::max(temp -> vTwo[1], std::max(temp -> vThree[1], maxY)));
                        Faces[chunkIndex].push_back(temp);
                    }
                }
            }
        }
        
        mainBoundary -> minimum.set(1, std::min(mainBoundary -> minimum[1], minY + cornerOffset[1]));
        mainBoundary -> maximum.set(1, std::max(mainBoundary -> maximum[1], maxY + cornerOffset[1]));
        boundaries.push_back(new BoundingBox(Vector(chunkCoord[0], minY, chunkCoord[2]) + cornerOffset, 
                                             Vector(chunkCoord[0] + cubesInChunk[0], maxY, chunkCoord[2] + cubesInChunk[2]) + cornerOffset));
    }

    public:
    Terrain() {}
    Terrain(const Vector& cornerPos, const Vector &tDim, const int &height, SimplexNoise* noise, const Vector& camOrigin, const short cPC = 32) {
        cornerOffset = cornerPos;
        Dimensions = tDim;
        maxHeight = height;
        halfHeight = tDim[1] * 0.5;
        cubesPerChunk = cPC;
        noiseMap = noise;
        Faces = std::vector<std::vector<Triangle*>>(0);
        boundaries = std::vector<BoundingBox*>(0);
        mainBoundary = new BoundingBox(Vector(0, INFINITY, 0), Vector(tDim[0], -1E16, tDim[2]));
        Generate(camOrigin);
    }

    void Generate(const Vector& camOrigin) {
        short xLimit, yLimit, zLimit;
        int chunkIndex = 0;
        for (short x = 0; x < Dimensions[0] - 1; x += cubesPerChunk) {
            xLimit = (Dimensions[0] - 1) - x >= cubesPerChunk ? cubesPerChunk : Dimensions[0] - 1 - x;
            for (short y = 0; y < Dimensions[1] - 1; y += cubesPerChunk) {
                yLimit = (Dimensions[1] - 1) - y >= cubesPerChunk ? cubesPerChunk : Dimensions[1] - 1 - y;
                for (short z = 0; z < Dimensions[2] - 1; z += cubesPerChunk) {
                    zLimit = (Dimensions[2] - 1) - z >= cubesPerChunk ? cubesPerChunk : Dimensions[2] - 1 - z;
                    Faces.push_back(std::vector<Triangle*>(0));
                    GenerateChunk(Vector(x, y, z), Vector(xLimit, yLimit, zLimit), chunkIndex++, camOrigin);
                }
            }
        }


        int numTriangles = 0;
        for (int i = 0; i < Faces.size(); ++i) {
            numTriangles += Faces[i].size();
        }

        std::cout << "Main Boundary: " << mainBoundary -> minimum << " <= P <= " << mainBoundary -> maximum << '\n';
        std::cout << "# of Boundaries: " << boundaries.size() << "\n";
        std::cout << "# of Triangles: " << numTriangles << "\n";
        std::cout << "Terrain Generated!\n";
    }

    bool intersect(const Ray& r, const double& min, const double& max, Record& info) const {
        if (!mainBoundary -> intersect(r, min, max)) return false;

        bool hitObject = false;
        double closestDistance = max;
        Record temp;

        for (int i = 0; i < boundaries.size(); ++i) {
            if (boundaries[i] -> intersect(r, min, closestDistance)) {
                for (Triangle* tri : Faces[i]) {
                    if (tri -> intersect(r, min, closestDistance, temp)) {
                        closestDistance = temp.distance;
                        hitObject = true;
                        info = temp;
                    }
                }
            }
        }

        return hitObject;
    }

    void buildObject(const Types& json) {
        Types* temp = json.jsonValue.get("dimensions");
        if (temp != nullptr) 
            Dimensions = *(temp -> vecValue);
        
        temp = json.jsonValue.get("corner_position");
        if (temp != nullptr)
            cornerOffset = *(temp -> vecValue);

        temp = json.jsonValue.get("maximum_height");
        if (temp != nullptr) 
            maxHeight = temp -> dblValue;
        
        temp = json.jsonValue.get("cubes_per_chunk");
        if (temp != nullptr)
            cubesPerChunk = (short) temp -> dblValue;
        else    
            cubesPerChunk = 32;

        halfHeight = Dimensions[1] * 0.5;
        Faces = std::vector<std::vector<Triangle*>>(0);
        boundaries = std::vector<BoundingBox*>(0);
        mainBoundary = new BoundingBox(Vector(0, 1E8, 0), Vector(Dimensions[0], -1E8, Dimensions[2]));
        material = buildMaterial(json);
        noiseMap = build3DNoise(json);
    }

    void printType() const {
        std::cout << "Terrain\n";
    }

};

#endif