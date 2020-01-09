#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <utility>

template <std::size_t x, std::size_t y>
using matrix = std::array<std::array<double, y>, x>;

constexpr int nx = 400;
constexpr int ny = 90;
constexpr int i1 = 200;
constexpr int i2 = 210;
constexpr int my_j1 = 50;
constexpr double delta = 0.01;
constexpr double sigma = 10 * delta;
constexpr double xA = 0.45;
constexpr double yA = 0.45;

constexpr double x(int i) {
    return delta * i;
}

constexpr double y(int j) {
    return delta * j;
}

constexpr matrix<nx + 1, ny + 1> fill_u0() {
    matrix<nx + 1, ny + 1> result = { { 0 } };
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            result[i][j] = 1.0 / (2 * M_PI * pow(sigma, 2)) * exp(-(pow(x(i) - xA, 2) + pow(y(j) - yA, 2)) / (2 * pow(sigma, 2)));
        }
    }
    return result;
}

constexpr std::pair<matrix<nx + 1, ny + 1>, matrix<nx + 1, ny + 1>> fill_v(const matrix<nx + 1, ny + 1>& psi) {
    matrix<nx + 1, ny + 1> vx = { { 0 } };
    matrix<nx + 1, ny + 1> vy = { { 0 } };
    for (int i = 1; i <= nx - 1; i++) {
        for (int j = 1; j <= ny - 1; j++) {
            vx[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2.0 * delta);
            vy[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2.0 * delta);
        }
    }
    for (int i = i1; i <= i2; i++) {
        for (int j = 0; j <= my_j1; j++) {
            vx[i][j] = 0.0;
            vy[i][j] = 0.0;
        }
    }
    for (int i = 1; i <= nx - 1; i++) {
        vx[i][0] = 0;
        vx[i][ny] = 0;
        vy[i][0] = 0;
        vy[i][ny] = 0;
    }
    for (int j = 0; j <= ny; j++) {
        vx[0][j] = vx[1][j];
        vx[nx][j] = vx[nx - 1][j];
    }

    return std::pair(vx, vy);
}

matrix<nx + 1, ny + 1> fill_psi() {
    matrix<nx + 1, ny + 1> result = { { 0 } };
    FILE* f = fopen("psi.dat", "r");
    int x = 0;
    int y = 0;
    double val = 0.0;
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            fscanf(f, "%d %d %lf", &x, &y, &val);
            result[x][y] = val;
        }
    }

    fclose(f);
    return result;
}

void iter_p(matrix<nx + 1, ny + 1>& u0, matrix<nx + 1, ny + 1>& u1, const matrix<nx + 1, ny + 1>& vx, const matrix<nx + 1, ny + 1>& vy, const double dt, const double D, const int it_max) {
    FILE* f = fopen(("out_" + std::to_string(D) + " " + std::to_string(dt) + ".dat").c_str(), "w");
    for (int it = 1; it <= it_max; it++) {
        for (int i = 0; i < u0.size(); i++) {
            std::memcpy(u1[i].data(), u0[i].data(), u0[i].size() * sizeof(double));
        }

        for (int k = 1; k <= 20; k++) {
            for (int i = 0; i <= nx; i++) {
                for (int j = 1; j <= ny - 1; j++) {
                    if (i > i1 && i < i2 && j < my_j1) {
                        continue;
                    }
                    if (i == 0) {
                        u1[i][j] = (1.0 / (1 + (2 * D * dt / pow(delta, 2)))) * (u0[i][j] - (dt / 2.0) * vx[i][j] * (((u0[i + 1][j] - u0[nx][j]) / (2.0 * delta)) + (u1[i + 1][j] - u1[nx][j]) / (2.0 * delta)) - (dt / 2) * vy[i][j] * ((u0[i][j + 1] - u0[i][j - 1]) / (2.0 * delta) + (u1[i][j + 1] - u1[i][j - 1]) / (2.0 * delta)) + dt / 2.0 * D * ((u0[i + 1][j] + u0[nx][j] + u0[i][j + 1] + u0[i][j - 1] - 4 * u0[i][j]) / pow(delta, 2) + (u1[i + 1][j] + u1[nx][j] + u1[i][j + 1] + u1[i][j - 1]) / pow(delta, 2)));
                    } else if (i == nx) {
                        u1[i][j] = (1.0 / (1 + (2 * D * dt / pow(delta, 2)))) * (u0[i][j] - (dt / 2.0) * vx[i][j] * (((u0[0][j] - u0[i - 1][j]) / (2.0 * delta)) + (u1[0][j] - u1[i - 1][j]) / (2.0 * delta)) - (dt / 2) * vy[i][j] * ((u0[i][j + 1] - u0[i][j - 1]) / (2.0 * delta) + (u1[i][j + 1] - u1[i][j - 1]) / (2.0 * delta)) + dt / 2.0 * D * ((u0[0][j] + u0[i - 1][j] + u0[i][j + 1] + u0[i][j - 1] - 4 * u0[i][j]) / pow(delta, 2) + (u1[0][j] + u1[i - 1][j] + u1[i][j + 1] + u1[i][j - 1]) / pow(delta, 2)));
                    } else {
                        u1[i][j] = (1.0 / (1 + (2 * D * dt / pow(delta, 2)))) * (u0[i][j] - (dt / 2.0) * vx[i][j] * (((u0[i + 1][j] - u0[i - 1][j]) / (2.0 * delta)) + (u1[i + 1][j] - u1[i - 1][j]) / (2.0 * delta)) - (dt / 2) * vy[i][j] * ((u0[i][j + 1] - u0[i][j - 1]) / (2.0 * delta) + (u1[i][j + 1] - u1[i][j - 1]) / (2.0 * delta)) + dt / 2.0 * D * ((u0[i + 1][j] + u0[i - 1][j] + u0[i][j + 1] + u0[i][j - 1] - 4 * u0[i][j]) / pow(delta, 2) + (u1[i + 1][j] + u1[i - 1][j] + u1[i][j + 1] + u1[i][j - 1]) / pow(delta, 2)));
                    }
                }
            }
        }
        for (int i = 0; i < u0.size(); i++) {
            std::memcpy(u0[i].data(), u1[i].data(), u0[i].size() * sizeof(double));
        }
        double c = 0.;
        double xsr = 0.;
        for (int i = 0; i <= nx; i++) {
            for (int j = 0; j <= ny; j++) {
                c += u0[i][j];
                xsr += x(i) * u0[i][j];
            }
        }
        c *= (delta * delta);
        xsr *= (delta * delta);
        fprintf(f, "%f %lf %lf\n", it * dt, c, xsr);
    }

    fclose(f);
}

void solution(double D, const int it_max) {
    auto u0 = fill_u0();
    matrix<nx + 1, ny + 1> u1 = { { 0 } };
    auto psi = fill_psi();
    auto [vx, vy] = fill_v(psi);
    double v_max = 0;
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            if (std::pow(vx[i][j], 2) + pow(vy[i][j], 2) > v_max) {
                v_max = std::pow(vx[i][j], 2) + pow(vy[i][j], 2);
            }
        }
    }
    v_max = std::sqrt(v_max);
    const double dt = delta / (4 * v_max);
    iter_p(u0, u1, vx, vy, dt, D, it_max);
}

int main() {
    solution(0, 1000);
    solution(0.1, 1000);
    return 0;
}