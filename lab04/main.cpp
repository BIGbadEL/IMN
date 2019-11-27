#include <array>
#include <cmath>
#include <functional>
#include <iostream>

constexpr double epsilon = 1.0;
constexpr double delta = 0.1;
constexpr std::size_t nx = 150;
constexpr std::size_t ny = 100;
constexpr double xmax = delta * nx;
constexpr double ymax = delta * ny;
constexpr double dx = 0.1 * xmax;
constexpr double dy = 0.1 * ymax;

template <std::size_t a, std::size_t b>
using my_array = std::array<std::array<double, b>, a>;

my_array<nx + 1, ny + 1> fill_array(double v1, double val, double v2) {
    my_array<nx + 1, ny + 1> result;
    for (auto& array : result) {
        for (auto& element : array) {
            element = val;
        }
    }
    for (int i = 0; i < nx + 1; i++) {
        result[i][0] = v1;
        result[i][ny] = v2;
    }
    return result;
}

constexpr double p1_val(double x, double y) {
    return std::exp(-std::pow(x - 0.35 * xmax, 2) / (dx * dx) - std::pow(y - 0.5 * ymax, 2) / (dy * dy));
}

constexpr double p2_val(double x, double y) {
    return -std::exp(-std::pow(x - 0.65 * xmax, 2) / (dx * dx) - std::pow(y - 0.5 * ymax, 2) / (dy * dy));
}

constexpr double p_val(double x, double y) {
    return p1_val(x, y) + p2_val(x, y);
}

my_array<nx, ny> fillP() {
    my_array<nx, ny> result;
    for (std::size_t i = 0; i < nx; i++) {
        for (std::size_t j = 0; j < ny; j++) {
            result[i][j] = p_val(i * delta, j * delta);
        }
    }
    return result;
}

constexpr double nextV(const my_array<nx, ny>& p,
    const my_array<nx + 1, ny + 1>& v, std::size_t i,
    std::size_t j) {
    return 0.25 * (v[i + 1][j] + v[i - 1][j] + v[i][j + 1] + v[i][j - 1] + delta * delta * p[i][j] / epsilon);
}

my_array<nx + 1, ny + 1> fill_new_V(my_array<nx + 1, ny + 1>& act_v,
    my_array<nx + 1, ny + 1>& old_v,
    const my_array<nx, ny>& p, double wg) {
    my_array<nx + 1, ny + 1> result = act_v;
    for (std::size_t i = 1; i < nx; i++) {
        for (std::size_t j = 1; j < ny; j++) {
            result[i][j] = nextV(p, old_v, i, j);
        }
    }
    for (std::size_t j = 0; j < ny + 1; j++) {
        result[0][j] = result[1][j];
        result[nx][j] = result[nx - 1][j];
    }
    for (std::size_t i = 0; i < nx + 1; i++) {
        for (std::size_t j = 0; j < ny + 1; j++) {
            old_v[i][j] = (1 - wg) * old_v[i][j] + wg * result[i][j];
        }
    }
    return result;
}

constexpr double calculate_S(const my_array<nx, ny>& p,
    const my_array<nx + 1, ny + 1>& v_n) {
    double result = 0.0;
    for (std::size_t i = 0; i < nx; i++) {
        for (std::size_t j = 0; j < ny; j++) {
            result += (delta * delta) * (0.5 * std::pow((v_n[i + 1][j] - v_n[i][j]) / delta, 2) + 0.5 * std::pow((v_n[i][j + 1] - v_n[i][j]) / delta, 2) - p[i][j] * v_n[i][j]);
        }
    }
    return result;
}

void global_relaxation(const my_array<nx, ny>& p, double wg, double TOL) {
    my_array<nx + 1, ny + 1> next_v = fill_array(10.0, 0.0, 0.0);
    my_array<nx + 1, ny + 1> start_v = fill_array(10.0, 0.0, 0.0);
    double Sit_prev = 0.0;
    double Sit = calculate_S(p, start_v);
    FILE* file = fopen(("1_" + std::to_string(wg) + ".dat").c_str(), "w");
    int i = 0;
    do {
        Sit_prev = Sit;
        next_v = fill_new_V(next_v, start_v, p, wg);
        Sit = calculate_S(p, next_v);
        fprintf(file, "%d %f\n", i++, Sit);
    } while (std::abs((Sit - Sit_prev) / Sit_prev) > TOL);
    fclose(file);
    file = fopen(("error_1_" + std::to_string(wg) + ".dat").c_str(), "w");
    for (int i = 1; i < nx; i++) {
        for (int j = 1; j < ny; j++) {
            double e = ((next_v[i + 1][j] - 2 * next_v[i][j] + next_v[i - 1][j]) / (delta * delta) + (next_v[i][j + 1] - 2 * next_v[i][j] + next_v[i][j - 1]) / (delta * delta)) + p[i][j] / epsilon;
            fprintf(file, "%f %f %f\n", i * delta, j * delta, e);
        }
    }
    fclose(file);
}

double nextV_local(double wl, const my_array<nx, ny>& p,
    const my_array<nx + 1, ny + 1>& v, int i, int j) {
    return (1 - wl) * v[i][j] + wl / 4.0 * (v[i + 1][j] + v[i - 1][j] + v[i][j + 1] + v[i][j - 1] + delta * delta / epsilon * p[i][j]);
}

void local_relaxation(const my_array<nx, ny>& p, double wl, double TOL) {
    my_array<nx + 1, ny + 1> v = fill_array(10.0, 0.0, 0.0);
    double Sit_prev = 0.0;
    double Sit = calculate_S(p, v);
    FILE* file = fopen(("2_" + std::to_string(wl) + ".dat").c_str(), "w");
    int i = 0;
    do {
        Sit_prev = Sit;
        for (int i = 1; i < nx; i++) {
            for (int j = 1; j < ny; j++) {
                v[i][j] = nextV_local(wl, p, v, i, j);
            }
        }
        for (int j = 1; j < ny; j++) {
            v[0][j] = v[1][j];
            v[nx][j] = v[nx - 1][j];
        }
        Sit = calculate_S(p, v);
        fprintf(file, "%d %f\n", i++, Sit);
    } while (std::abs((Sit - Sit_prev) / Sit_prev) > TOL);
    fclose(file);
    file = fopen(("error_2_" + std::to_string(wl) + ".dat").c_str(), "w");
    for (int i = 1; i < nx; i++) {
        for (int j = 1; j < ny; j++) {
            double e = ((v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) / (delta * delta) + (v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]) / (delta * delta)) + p[i][j] / epsilon;
            fprintf(file, "%f %f %f\n", i * delta, j * delta, e);
        }
    }
    fclose(file);
}

auto main() -> int {
    my_array<nx, ny> p = fillP();
    global_relaxation(p, 0.6, 1e-8);
    global_relaxation(p, 1.0, 1e-8);
    local_relaxation(p, 1.0, 1e-8);
    local_relaxation(p, 1.4, 1e-8);
    local_relaxation(p, 1.8, 1e-8);
    local_relaxation(p, 1.9, 1e-8);
    return 0;
}

// template <typename T, std::size_t size>
// class arr {
// public:
//     constexpr arr() {
//         for(int i = 0; i < size; i++){
//             m_data[i] = T{};
//         }
//     }
//     T& operator[] (std::size_t index) {
//         return m_data[index];
//     }
//     T operator[] (std::size_t index) const {
//         return m_data[index];
//     }
// private:
//     std::size_t m_size = size;
//     T m_data[size];
// };

// constexpr arr<arr<double, nx>, ny> create() {
//     arr<arr<double, nx>, ny> result;
//     for(int i = 0; i < nx; i++){
//         for(int j = 0; j < ny; j++){
//             result[i][j] = 1.0;
//         }
//     }
//     return result;
// }