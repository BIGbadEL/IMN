#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

struct Point {
    double x;
    double y;
};
constexpr int p = 2;
constexpr double x0 = 0.01;
constexpr double v0 = 0.0;
constexpr double dt0 = 1.0;
constexpr double S = 0.75;
constexpr int tmax = 40;
constexpr double alfa = 5.0;

constexpr double f(double v) { return v; }
constexpr double g(double x, double v) { return alfa * (1 - x * x) * v - x; }

Point trap_next(double xn, double vn, double dt, double alfa_local) {
    constexpr double a11 = 1;
    double dx = 0.0, dv = 0.0;
    double x_next = xn;
    double v_next = vn;
    do {
        const double a12 = -dt / 2;
        const double a21 = -dt / 2.0 * (-2.0 * alfa_local * x_next * v_next - 1);
        const double a22 = 1 - dt / 2 * alfa_local * (1 - x_next * x_next);
        const double F = x_next - xn - dt / 2 * (f(vn) + f(v_next));
        const double G = v_next - vn - dt / 2 * (g(xn, vn) + g(x_next, v_next));
        dx = ((-F) * a22 - (-G) * a12) / (a11 * a22 - a12 * a21);
        dv = ((-G) * a11 - (-F) * a21) / (a11 * a22 - a12 * a21);
        x_next += dx;
        v_next += dv;
    } while (std::abs(dx) > 1e-10 && std::abs(dv) > 1e-10);
    return Point { x_next, v_next };
}

constexpr Point RK2_next(double xn, double vn, double dt, double alfa_local) {
    const double k1x = vn;
    const double k1v = alfa_local * (1 - xn * xn) * vn - xn;
    const double k2x = vn + dt * k1v;
    const double k2v = g(xn + dt * k1x, vn + dt * k1v);
    return Point { xn + dt / 2 * (k1x + k2x), vn + dt / 2 * (k1v + k2v) };
}

constexpr double new_dt(double old_dt, double Ex, double Ev, double tol) {
    double E_max = std::max(Ex, Ev);
    return std::pow((S * tol) / E_max, 1.0 / (p + 1.0)) * old_dt;
}

constexpr double error(double x_2, double x_1, double loc_p) {
    return std::abs((x_2 - x_1) / (std::pow(2, loc_p) - 1));
}

template <typename callable>
void solution(double tol, const char* name, callable fun) {
    double t = 0.0;
    double dt = dt0;
    double xn = x0;
    double vn = v0;
    FILE* file = fopen(name, "wb");
    do {
        Point p_2 = fun(xn, vn, dt, alfa);
        p_2 = fun(p_2.x, p_2.y, dt, alfa);

        Point p_1 = fun(xn, vn, 2 * dt, alfa);
        double Ex = error(p_2.x, p_1.x, p);
        double Ev = error(p_2.y, p_1.y, p);
        if (std::max(Ex, Ev) < tol) {
            t += 2 * dt;
            xn = p_2.x;
            vn = p_2.y;
            fprintf(file, "%f %f %f %f\n", t, dt, p_2.x, p_2.y);
        }
        dt = new_dt(dt, Ex, Ev, tol);
    } while (t < tmax);
    fclose(file);
}

int main() {
    solution(1e-2, "zad1.dat", trap_next);
    solution(1e-5, "zad1b.dat", trap_next);
    solution(1e-2, "zad2.dat", RK2_next);
    solution(1e-5, "zad2b.dat", RK2_next);
    return 0;
}