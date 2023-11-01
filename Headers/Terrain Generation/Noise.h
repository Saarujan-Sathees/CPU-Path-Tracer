#ifndef NOISE
#define NOISE

#include "../Vector.h"
#include "../Random.h"
#include "../JSON.h"
#include "LookupTable.h"
#include <vector>
#include <iostream>
#include <future>
#include <atomic>
#include <mutex>

using Arr = std::vector<double>;
using Arr2 = std::vector<Arr>;
using Arr3 = std::vector<Arr2>;

class SimplexNoise {
    double nSum, s, t0, t1, t2, t3, amplitude;
    char gi0, gi1, gi2, gi3;
    int nPT;
    Vector zero, one, two, three;
    int_fast8_t i, j, k, ii, jj, kk, i1, j1, k1, i2, j2, k2, i3, j3, k3;

    const double G3 = 1.0 / 6.0, G2 = G3 * 2.0, G1 = G3 * 3.0;
    Arr3 values;
    Vector worldLimit;
    float minMaxAv, negMinMax;

    inline int fastFloor(const double &d) const {
        return d > 0 ? (int) d : int(d - 1);
    }

    public:
    SimplexNoise(const double &amp, const int &mi, const int &ma, const int noisePerThread = 0) {
        amplitude = amp;
        minMaxAv = (mi + ma) / 2.0;
        negMinMax = (ma - mi) / 2.0;
        nPT = noisePerThread;
    };

    inline void Generate(const Vector &wL, const double &scale, const double &persistence, const double &lacunarity, const int &iterations = 3) {     
        const std::clock_t begin = std::clock();
        worldLimit = wL;

        if (nPT > 0) {
            std::atomic<int> Count(0);
            std::vector<std::future<void>> Threads;
            std::mutex m;
            values = Arr3(worldLimit[0], Arr2(worldLimit[1], Arr(worldLimit[2])));

            std::cout << "Values Made\n";
            for (int Thread = 0; Thread < worldLimit[0]; Thread += nPT) {
                Threads.emplace_back(std::async(std::launch::deferred, [=, &Count, &m]() {
                    const int xStart = Count;
                    Count += nPT;
                    Arr3 xyz((wL[0] - xStart < nPT ? wL[0] - xStart : nPT), Arr2(worldLimit[1], Arr(worldLimit[2])));
                    
                    for (int x = 0; x < xyz.size(); ++x) {
                        for (int y = 0; y < xyz[x].size(); ++y) {
                            for (int z = 0; z < xyz[x][y].size(); ++z) {
                                xyz[x][y][z] = noiseSum(Vector(x + xStart, y, z), scale, persistence, lacunarity, iterations);
                            }
                        }
                    }
                    
                    m.lock();
                    for (int i = 0; i < xyz.size(); ++i) {
                        values[xStart + i] = xyz[i];
                    }
                    m.unlock();
                })); 
                
            }

            for (int i = 0; i < Threads.size(); ++i) Threads[i].wait();
            m.unlock();
            
        } else {
            values = Arr3(worldLimit[0], Arr2(worldLimit[1], Arr(worldLimit[2])));
            for (int x = 0; x < worldLimit[0]; ++x) {
                for (int y = 0; y < worldLimit[1]; ++y) {
                    for (int z = 0; z < worldLimit[2]; ++z) {
                        values[x][y][z] = noiseSum(Vector(x, y, z), scale, persistence, lacunarity, iterations);
                    }
                }
            }
        }

        
        std::cout << "Noise Time: " << int(double(std::clock() - begin) / CLOCKS_PER_SEC) << "\n";
    }

    inline double valueAt(const Vector &xyz) const {
        return values[xyz[0]][xyz[1]][xyz[2]];
    }

    inline double noiseSum(const Vector &xyz, const double &scale, const double &persistence, const double &lacunarity, const int &iterations = 3) {
        double maxAmplitude = 0, noise = 0, freq = scale, amp = amplitude;

        for (char i = 0; i < iterations; ++i) {
            noise += Noise(xyz * freq) * amp;
            maxAmplitude += amp;
            amp *= persistence;
            freq *= lacunarity;
        }

        //noise /= maxAmplitude;
        return noise * negMinMax + minMaxAv;
    }

    inline double Noise(const Vector &xyz) {
        nSum = 0;
        s = (xyz[0] + xyz[1] + xyz[2]) * G2; //Skew Factor
        i = fastFloor(xyz[0] + s), j = fastFloor(xyz[1] + s), k = fastFloor(xyz[2] + s);
        zero = Vector(xyz[0] - i, xyz[1] - j, xyz[2] - k) + (i + j + k) * G3;

        if (zero[0] >= zero[1]) {
            j1 = 0;
            i2 = 1;
            if (zero[1] >= zero[2]) { //X Y Z
                i1 = 1;  
                k1 = 0; 
                j2 = 1; 
                k2 = 0; 
            } else if (zero[0] >= zero[2]) { //X Z Y
                i1 = 1;  
                k1 = 0;  
                j2 = 0; 
                k2 = 1; 
            } else { //Z X Y
                i1 = 0;  
                k1 = 1; 
                j2 = 0; 
                k2 = 1; 
            } 
        } else { 
            i1 = 0;
            j2 = 1;
            if (zero[1] < zero[2]) { //Z Y X
                j1 = 0; 
                k1 = 1; 
                i2 = 0;  
                k2 = 1; 
            } else if (zero[0] < zero[2]) { //Y Z X
                j1 = 1; 
                k1 = 0; 
                i2 = 0; 
                k2 = 1; 
            } else { //Y X Z
                j1 = 1; 
                k1 = 0; 
                i2 = 1; 
                k2 = 0; 
            } 
        }
    
        one = Vector(zero[0] - i1, zero[1] - j1, zero[2] - k1) + G3;
        two = Vector(zero[0] - i2, zero[1] - j2, zero[2] - k2) + G2;
        three = Vector(zero[0], zero[1], zero[2]) - 1 + G1;
        
        ii = i & 255, jj = j & 255, kk = k & 255;

        gi0 = perm[ii + perm[jj + perm[kk]]] % 12;
        gi1 = perm[ii + i1 + perm[jj + j1 + perm[kk + k1]]] % 12;
        gi2 = perm[ii + i2 + perm[jj + j2 + perm[kk + k2]]] % 12;
        gi3 = perm[ii + 1 + perm[jj + 1 + perm[kk + 1]]] % 12;
        
        t0 = 0.6 - zero.squaredLength();
        if (t0 >= 0) 
            nSum += std::pow(t0, 4) * zero.dot(gradient[gi0]);

        t1 = 0.6 - one.squaredLength();
        if (t1 >= 0) 
            nSum += std::pow(t1, 4) * one.dot(gradient[gi1]);
        
        t2 = 0.6 - two.squaredLength();
        if (t2 >= 0)
            nSum += std::pow(t2, 4) * two.dot(gradient[gi2]);

        t3 = 0.6 - three.squaredLength();
        if (t3 >= 0) 
            nSum += std::pow(t3, 4) * three.dot(gradient[gi3]);
        
        return 32.0 * nSum;
    }

};

class Noise2D {
    double amplitude;
    int nPT, length, width;
    std::vector<Vector> octaveOffsets;
    const double randMult = (3.14159265 / ~(~0u >> 1));

    Arr2 values;
    float minMaxAv, negMinMax;

    inline int fastFloor(const double &d) const {
        return d > 0 ? (int) d : int(d - 1);
    }

    public:
    Noise2D(const double &amp, const int &mi, const int &ma, const int noisePerThread = 0) {
        amplitude = amp;
        nPT = noisePerThread;
    };

    inline void Generate(const Vector &wL, const double &scale, const double &pers, const double &lac, const int &iter = 3, const int& seed = (int) randDbl(0, 1E6)) {     
        const std::clock_t begin = std::clock();
        length = (int) wL[0];
        width = (int) wL[2];

        std::minstd_rand reng(seed);
        std::uniform_real_distribution<double> randDist(wL[0], wL[0] + 1);
        FILE* output = fopen("noise.ppm", "w");
        fprintf(output, "P3\n%d %d\n255\n", length, width);

        for (int i = 0; i < iter; ++i) {
            octaveOffsets.push_back(Vector((int) randDist(reng), 0, (int) randDist(reng)));
        }

        if (nPT > 0) {
            std::atomic<int> Count(0);
            std::atomic<double> min(INFINITY), max(-1E8);
            std::vector<std::future<void>> Threads;
            std::mutex m;
            values = Arr2(length, Arr(width));

            for (int Thread = 0; Thread < length; Thread += nPT) {
                Threads.emplace_back(std::async(std::launch::deferred, [=, &Count, &m, &min, &max]() {
                    const int xStart = Count;
                    Count += nPT;
                    Arr2 xz(std::min(length - xStart, nPT), Arr(width));
                    
                    for (int x = 0; x < xz.size(); ++x) {
                        for (int z = 0; z < xz[x].size(); ++z) {
                            xz[x][z] = noiseSum(x + xStart, z, scale, pers, lac, iter);
                            if (xz[x][z] < min)
                                min = values[x][z];
                            
                            if (values[x][z] > max)
                                max = values[x][z];
                        }
                    }
                    
                    m.lock();
                    for (int i = 0; i < xz.size(); ++i) {
                        values[xStart + i] = xz[i];
                    }
                    m.unlock();
                })); 
                
            }

            for (int i = 0; i < Threads.size(); ++i) Threads[i].wait();
            max = max - min;

            for (int x = 0; x < length; ++x) {
                for (int z = 0; z < width; ++z) {
                    double clamped = (values[x][z] - min) / max * 255;
                    fprintf(output, "%lf %lf %lf\n", clamped, clamped, clamped);
                }
            }

            fclose(output);
            
        } else {
            double min = INFINITY, max = -1E8;
            values = Arr2(length, Arr(width));
            for (int x = 0; x < length; ++x) {
                for (int z = 0; z < width; ++z) {
                    values[x][z] = noiseSum(x, z, scale, pers, lac, iter);
                    if (values[x][z] < min)
                        min = values[x][z];
                    
                    if (values[x][z] > max)
                        max = values[x][z];
                }
            }

            max -= min;

            for (int x = 0; x < length; ++x) {
                for (int z = 0; z < width; ++z) {
                    double clamped = (values[x][z] - min) / max * 255;
                    fprintf(output, "%lf %lf %lf\n", clamped, clamped, clamped);
                }
            }

            fclose(output);
        }

        std::cout << "Noise Time: " << int(double(std::clock() - begin) / CLOCKS_PER_SEC) << "\n";
    }

    inline double valueAt(const double& x, const double& z) const {
        return values[x][z];
    }

    inline double noiseSum(const double& x, const double& z, const double &scale, const double &pers, const double &lac, const int &iter = 3) {
        double noise = 0, freq = scale, amp = amplitude;
        Vector position(x, 0, z);

        for (int i = 0; i < iter; ++i) {
            noise += generateNoise(x * freq + octaveOffsets[i][0], z * freq + octaveOffsets[i][2]) * amp;
            amp *= pers;
            freq *= lac;
        }

        //noise /= maxAmplitude;
        return noise;
    }

    /*
    inline double generateNoise(const Vector &xyz) {
        nSum = 0;
        s = (xyz[0] + xyz[1] + xyz[2]) * G2; //Skew Factor
        i = fastFloor(xyz[0] + s), j = fastFloor(xyz[1] + s), k = fastFloor(xyz[2] + s);
        zero = Vector(xyz[0] - i, xyz[1] - j, xyz[2] - k) + (i + j + k) * G3;

        if (zero[0] >= zero[1]) {
            j1 = 0;
            i2 = 1;
            if (zero[1] >= zero[2]) { //X Y Z
                i1 = 1;  
                k1 = 0; 
                j2 = 1; 
                k2 = 0; 
            } else if (zero[0] >= zero[2]) { //X Z Y
                i1 = 1;  
                k1 = 0;  
                j2 = 0; 
                k2 = 1; 
            } else { //Z X Y
                i1 = 0;  
                k1 = 1; 
                j2 = 0; 
                k2 = 1; 
            } 
        } else { 
            i1 = 0;
            j2 = 1;
            if (zero[1] < zero[2]) { //Z Y X
                j1 = 0; 
                k1 = 1; 
                i2 = 0;  
                k2 = 1; 
            } else if (zero[0] < zero[2]) { //Y Z X
                j1 = 1; 
                k1 = 0; 
                i2 = 0; 
                k2 = 1; 
            } else { //Y X Z
                j1 = 1; 
                k1 = 0; 
                i2 = 1; 
                k2 = 0; 
            } 
        }
    
        one = Vector(zero[0] - i1, zero[1] - j1, zero[2] - k1) + G3;
        two = Vector(zero[0] - i2, zero[1] - j2, zero[2] - k2) + G2;
        three = Vector(zero[0], zero[1], zero[2]) - 1 + G1;
        
        ii = i & 255, jj = j & 255, kk = k & 255;

        gi0 = perm[ii + perm[jj + perm[kk]]] % 12;
        gi1 = perm[ii + i1 + perm[jj + j1 + perm[kk + k1]]] % 12;
        gi2 = perm[ii + i2 + perm[jj + j2 + perm[kk + k2]]] % 12;
        gi3 = perm[ii + 1 + perm[jj + 1 + perm[kk + 1]]] % 12;
        
        t0 = 0.6 - zero.squaredLength();
        if (t0 >= 0) 
            nSum += std::pow(t0, 4) * zero.dot(gradient[gi0]);

        t1 = 0.6 - one.squaredLength();
        if (t1 >= 0) 
            nSum += std::pow(t1, 4) * one.dot(gradient[gi1]);
        
        t2 = 0.6 - two.squaredLength();
        if (t2 >= 0)
            nSum += std::pow(t2, 4) * two.dot(gradient[gi2]);

        t3 = 0.6 - three.squaredLength();
        if (t3 >= 0) 
            nSum += std::pow(t3, 4) * three.dot(gradient[gi3]);
        
        return 32.0 * nSum;
    }
    */

    inline double interpolate(const double& a, const double& b, const double& t) const {
        if (t < 0) return a;
        if (t > 1) return b;

        return (b - a) * (3 - t * 2) * std::pow(t, 2) + a;
    }

    inline Vector randomGradient(const int& x, const int& y) const {
        const unsigned w = 8 * sizeof(unsigned);
        const unsigned s = w / 2; 
        unsigned a = x, b = y;
        a *= 3284157443; 
        b ^= a << s | a >> w-s;
        b *= 1911520717; 
        a ^= b << s | b >> w-s;
        a *= 2048419325;

        double random = a * randMult;
        return Vector(std::cos(random), std::sin(random), 0);
    }

    inline double dotGradient(const int& xA, const int& yA, const double& x, const double& y) const {
        Vector gradient = randomGradient(xA, yA);

        double distX = x - xA, distY = y - yA;
        return gradient[0] * distX + gradient[1] * distY;
    }

    inline double generateNoise(const double& x, const double& y) {
        int xOne = (int) x, yOne = (int) y, xTwo = xOne + 1, yTwo = yOne + 1;

        double xWeight = x - xOne, yWeight = y - yOne;
        double ixOne = interpolate(dotGradient(xOne, yOne, x, y), dotGradient(xTwo, yOne, x, y), xWeight),
               ixTwo = interpolate(dotGradient(xOne, yTwo, x, y), dotGradient(xTwo, yTwo, x, y), xWeight);

        return (interpolate(ixOne, ixTwo, yWeight));
    }
};

SimplexNoise* build3DNoise(const Types& json) {
    double amplitude, lacunarity, scale, persistence;
    Vector dimensions;
    int min, max, noisePerThread = 0, iterations;

    Types* temp = json.jsonValue.get("noise_amplitude");
    if (temp != nullptr)
        amplitude = temp -> dblValue;

    temp = json.jsonValue.get("noise_min");
    if (temp != nullptr)
        min = (int) temp -> dblValue;

    temp = json.jsonValue.get("noise_max");
    if (temp != nullptr)
        max = (int) temp -> dblValue;

    temp = json.jsonValue.get("noise_per_thread");
    if (temp != nullptr)
        noisePerThread = (int) temp -> dblValue;

    temp = json.jsonValue.get("noise_lacunarity");
    if (temp != nullptr)
        lacunarity = temp -> dblValue;

    temp = json.jsonValue.get("noise_scale");
    if (temp != nullptr)
        scale = temp -> dblValue;
    
    temp = json.jsonValue.get("noise_persistence");
    if (temp != nullptr)
        persistence = temp -> dblValue;

    temp = json.jsonValue.get("noise_octaves");
    if (temp != nullptr)
        iterations = (int) temp -> dblValue;

    temp = json.jsonValue.get("dimensions");
    if (temp != nullptr)
        dimensions = *(temp -> vecValue);

    SimplexNoise* res = new SimplexNoise(amplitude, min, max, noisePerThread);
    res -> Generate(dimensions, scale, persistence, lacunarity, iterations);
    return res;
}

Noise2D* build2DNoise(const Types& json) {
    double amplitude, lacunarity, scale, persistence;
    Vector dimensions;
    int min, max, noisePerThread = 0, iterations, seed = (int) randDbl(0, 1E6);

    Types* temp = json.jsonValue.get("noise_amplitude");
    if (temp != nullptr)
        amplitude = temp -> dblValue;

    temp = json.jsonValue.get("noise_min");
    if (temp != nullptr)
        min = (int) temp -> dblValue;

    temp = json.jsonValue.get("noise_max");
    if (temp != nullptr)
        max = (int) temp -> dblValue;

    temp = json.jsonValue.get("noise_per_thread");
    if (temp != nullptr)
        noisePerThread = (int) temp -> dblValue;

    temp = json.jsonValue.get("noise_lacunarity");
    if (temp != nullptr)
        lacunarity = temp -> dblValue;

    temp = json.jsonValue.get("noise_scale");
    if (temp != nullptr)
        scale = temp -> dblValue;
    
    temp = json.jsonValue.get("noise_persistence");
    if (temp != nullptr)
        persistence = temp -> dblValue;

    temp = json.jsonValue.get("noise_octaves");
    if (temp != nullptr)
        iterations = (int) temp -> dblValue;

    temp = json.jsonValue.get("dimensions");
    if (temp != nullptr)
        dimensions = *(temp -> vecValue);

    temp = json.jsonValue.get("map_seed");
    if (temp != nullptr)
        seed = (int) temp -> dblValue;

    Noise2D* res = new Noise2D(amplitude, min, max, noisePerThread);
    res -> Generate(dimensions, scale, persistence, lacunarity, iterations, seed);
    return res;
}

#endif