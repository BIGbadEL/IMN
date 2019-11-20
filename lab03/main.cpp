#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

struct Point {
    double x;
    double y;
};

struct my_file {
public:
    my_file(const char* name)
        : _fptr(fopen(name, "wb")) {}
    ~my_file() { fclose(_fptr); }
    void lab03_print(double t, double dt, Point& p) {
        fprintf(_fptr, "%f %f %f %f\n", t, dt, p.x, p.y);
    }

private:
    FILE* _fptr = nullptr;
};

constexpr int p = 2;
constexpr double x0 = 0.01;
constexpr double v0 = 0.0;
constexpr double dt0 = 1.0;
constexpr double S = 0.75;
constexpr int tmax = 40;
constexpr double alfa = 5.0;

constexpr Point trap_next(double xn, double vn, double dt, double alfa) {
    return Point { 0, 0 };
}

constexpr Point RK2_next(double xn, double vn, double dt, double alfa) {
    return Point { 0, 0 };
}

constexpr double new_dt(double old_dt, double Ex, double Ev, double tol) {
    double E_max = std::max(Ex, Ev);
    return std::pow((S * tol) / E_max, 1.0 / (p + 1.0)) * old_dt;
}

constexpr double error(double x_2, double x_1, double loc_p) {
    return (x_2 - x_1) / (std::pow(2, loc_p) - 1);
}

template <typename callable>
void solution(double tol, callable fun) {
    double t = 0.0;
    double dt = dt0;
    double xn = x0;
    double vn = v0;
    do {
        Point p_2 = fun(xn, vn, dt, alfa);
        p_2 = fun(p_2.x, p_2.y, dt, alfa);

        Point p_1 = fun(xn, vn, 2 * dt, alfa);
        double Ex = error(p_2.x, p_1.x, p);
        double Ev = error(p_2.y, p_1.y);
        if (std::max(Ex, Ev) < tol) {
            t += 2 * dt;
        }
        dt = new_dt(dt, Ex, Ev, tol);
    } while (t < tmax);
}

int main() {

    return 0;
}