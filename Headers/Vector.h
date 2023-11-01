#ifndef VECTOR
#define VECTOR

#include <math.h>
#include <ostream>

class Vector {
    protected:
    double values[3];

    public:
    Vector() {
        values[0] = 0;
        values[1] = 0;
        values[2] = 0;
    }

    Vector(double x, double y, double z) {
        values[0] = x;
        values[1] = y;
        values[2] = z;
    }

    #pragma region Operators
    double operator [](const int& index) const {
        return values[index];
    }

    void set(const int& index, const double& value) {
        values[index] = value;
    }

    Vector operator -() const {
        return Vector(-values[0], -values[1], -values[2]);
    }

    inline Vector operator + (const Vector &b) const {
        return Vector(values[0] + b.values[0], values[1] + b.values[1], values[2] + b.values[2]);
    }

    inline Vector operator + (const double &b) const {
        return Vector(values[0] + b, values[1] + b, values[2] + b);
    }

    inline Vector operator - (const Vector &b) const {
        return Vector(values[0] - b.values[0], values[1] - b.values[1], values[2] - b.values[2]);
    }

    inline Vector operator - (const double &b) const {
        return Vector(values[0] - b, values[1] - b, values[2] - b);
    }

    inline Vector operator * (const Vector &b) const {
        return Vector(values[0] * b.values[0], values[1] * b.values[1], values[2] * b.values[2]);
    }

    inline Vector operator * (const double &b) const {
        return Vector(values[0] * b, values[1] * b, values[2] * b);
    }

    inline Vector operator / (const Vector &b) const {
        return Vector(values[0] / b.values[0], values[1] / b.values[1], values[2] / b.values[2]);
    }

    inline Vector operator / (const double &b) const {
        const double reciprocal = 1.0 / b;
        return Vector(values[0] * reciprocal, values[1] * reciprocal, values[2] * reciprocal);
    }

    Vector& operator += (const Vector &b) {
        values[0] += b.values[0];
        values[1] += b.values[1];
        values[2] += b.values[2];
        return *this;
    }

    Vector& operator += (const double &b) {
        values[0] += b;
        values[1] += b;
        values[2] += b;
        return *this;
    }

    Vector& operator -= (const Vector &b) {
        values[0] -= b.values[0];
        values[1] -= b.values[1];
        values[2] -= b.values[2];
        return *this;
    }

    Vector& operator -= (const double &b) {
        values[0] -= b;
        values[1] -= b;
        values[2] -= b;
        return *this;
    }

    Vector& operator *= (const Vector &b) {
        values[0] *= b.values[0];
        values[1] *= b.values[1];
        values[2] *= b.values[2];
        return *this;
    }

    Vector& operator *= (const double &b) {
        values[0] *= b;
        values[1] *= b;
        values[2] *= b;
        return *this;
    }

    Vector& operator /= (const Vector &b) {
        values[0] /= b.values[0];
        values[1] /= b.values[1];
        values[2] /= b.values[2];
        return *this;
    }

    Vector& operator /= (const double &b) {
        values[0] /= b;
        values[1] /= b;
        values[2] /= b;
        return *this;
    }
    #pragma endregion

    inline double dot(const Vector& b) const {
        return values[0] * b.values[0] + values[1] * b.values[1] + values[2] * b.values[2];
    }

    inline Vector cross(const Vector& b) const {
        return Vector(values[1] * b.values[2] - b.values[1] * values[2],
                      -(values[0] * b.values[2] - b.values[0] * values[2]),
                      values[0] * b.values[1] - b.values[0] * values[1]);
    }

    inline Vector sqrt() const {
        return Vector(std::sqrt(values[0]), std::sqrt(values[1]), std::sqrt(values[2]));
    }

    inline Vector abs() const {
        return Vector(std::abs(values[0]), std::abs(values[1]), std::abs(values[2]));
    }

    inline double squaredLength() const {
        return values[0] * values[0] + values[1] * values[1] + values[2] * values[2];
    }

    inline double length() const {
        return std::sqrt(squaredLength());
    }

    inline Vector unitVector() const {
        return *this / length();
    }

    inline friend std::ostream& operator << (std::ostream& output, Vector v) {
        output << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
        return output;
    }

    bool nearZero() {
        static const double min = 1E-8;
        return std::fabs(values[0]) <= min && std::fabs(values[1]) <= min && std::fabs(values[2]) <= min; 
    }

    std::string str() {
        return std::to_string(values[0]) + " " + std::to_string(values[1]) + " " + std::to_string(values[2]);
    }

    inline Vector replaceNAN() const {
        return Vector(values[0] != values[0] ? 0 : values[0],
                      values[1] != values[1] ? 0 : values[1],
                      values[2] != values[2] ? 0 : values[2]);
    }
};

Vector interpolate(const Vector& a, const Vector& b, const double& amount) {
    return (b - a) * amount + a;//a * amount + b * (1.0 - amount);
} 

Vector clamp(const Vector& a, const double& min = 0, const double& max = 1) {
    return Vector(a[0] < min ? min : a[0] > max ? max : a[0],
                  a[1] < min ? min : a[1] > max ? max : a[1],
                  a[2] < min ? min : a[2] > max ? max : a[2]);
}

class Color : public Vector {
    const static double colorMult;
    public:
    Color() {
        values[0] = 0;
        values[1] = 0;
        values[2] = 0;
    }

    Color(const double& r, const double& g, const double& b) {
        values[0] = r * colorMult;
        values[1] = g * colorMult;
        values[2] = b * colorMult;
    }

    Color(const Vector& v) {
        values[0] = v[0] * colorMult;
        values[1] = v[1] * colorMult;
        values[2] = v[2] * colorMult;
    }
};

const Vector BLACK = Vector(0, 0, 0), WHITE = Color(255, 255, 255);
const double Color::colorMult = 1.0 / 255;

#endif