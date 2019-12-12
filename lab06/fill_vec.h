#pragma once
#include <cmath>
#include <iostream>
#include <vector>

constexpr double delta = 0.1;

double p_1(const double x, const double y, const double x_max, const double y_max) {
    return std::exp(-std::pow((x - 0.25 * x_max) / (0.1 * x_max), 2) - pow((y - 0.5 * y_max) / (0.1 * y_max), 2));
}
double p_2(const double x, const double y, const double x_max, const double y_max) {
    return -1.0 * exp(-pow((x - 0.75 * x_max) / (0.1 * x_max), 2) - pow((y - 0.5 * y_max) / (0.1 * y_max), 2));
}
double p(const double x, const double y, const double x_max, const double y_max) {
    return p_1(x, y, x_max, y_max) + p_2(x, y, x_max, y_max);
}

constexpr std::size_t calc_l(const std::size_t i, const std::size_t j, const std::size_t nx) {
    return i + j * (nx + 1);
}

constexpr std::size_t calc_j(const std::size_t l, const std::size_t nx) {
    return l / (nx + 1);
}

constexpr std::size_t calc_i(const std::size_t l, const std::size_t j, const std::size_t nx) {
    return l - j * (nx + 1);
}

void fill_vecs(std::size_t nx, std::size_t ny,
    const double e1, const double e2,
    std::vector<double>& b, std::vector<int>& ia,
    std::vector<double>& a, std::vector<int>& ja,
    int& nz_num, const double V1, const double V2,
    const double V3, const double V4, bool is_zero = true) {
    int k = -1;
    const std::size_t N = (nx + 1) * (ny + 1);

    b.clear();
    b.resize(N);

    ia.clear();
    ia.resize(N + 1);

    ja.clear();
    ja.resize(5 * N);

    a.clear();
    a.resize(5 * N);

    std::vector<double> e;
    for (std::size_t l = 0; l < N + nx + 1; l++) {
        const std::size_t i = calc_i(l, calc_j(l, nx), nx);
        e.push_back(i <= nx / 2 ? e1 : e2);
    }

    for (std::size_t l = 0; l < N; l++) {
        const std::size_t j = calc_j(l, nx);
        const std::size_t i = calc_i(l, j, nx);
        std::size_t brzeg = 0;
        double vb = 0.0;

        if (i == 0) {
            brzeg = 1;
            vb = V1;
        } else if (i == nx) {
            brzeg = 1;
            vb = V3;
        }

        if (j == ny) {
            brzeg = 1;
            vb = V2;
        } else if (j == 0) {
            brzeg = 1;
            vb = V4;
        }
        b[l] = (is_zero) ? 0 : -p(delta * i, delta * j, nx * delta, ny * delta);
        if (brzeg == 1) {
            b[l] = vb;
        }
        ia[l] = -1;
        if (l - nx - 1 >= 0 && brzeg == 0) {
            k++;
            if (ia[l] < 0) {
                ia[l] = k;
            }
            a[k] = e[l] / (delta * delta);
            ja[k] = l - nx - 1;
        }

        if (l - 1 >= 0 && brzeg == 0) {
            k++;
            if (ia[l] < 0) {
                ia[l] = k;
            }
            a[k] = e[l] / (delta * delta);
            ja[k] = l - 1;
        }

        k++;
        if (ia[l] < 0) {
            ia[l] = k;
        }

        if (brzeg == 0) {
            a[k] = -(2 * e[l] + e[l + 1] + e[l + nx + 1]) / (delta * delta); //fix
        } else {
            a[k] = 1;
        }
        ja[k] = l;
        if (l < N && brzeg == 0) {
            k++;
            a[k] = e[l + 1] / (delta * delta);
            ja[k] = l + 1;
        }

        if (l < N - nx - 1 && brzeg == 0) {
            k++;
            a[k] = e[l + nx + 1] / (delta * delta);
            ja[k] = l + nx + 1;
        }
    }
    nz_num = k + 1;
    ia[N] = nz_num;
    FILE* file = fopen(("a_" + std::to_string(nx) + "_" + std::to_string(e1) + ".dat").c_str(), "w");
    for (std::size_t l = 0; l < N; l++) {
        const std::size_t j = calc_j(l, nx);
        const std::size_t i = calc_i(l, j, nx);
        fprintf(file, "%ld %ld %ld %f\n", l, i, j, a[l]);
    }
    fclose(file);

    file = fopen(("b_" + std::to_string(nx) + "_" + std::to_string(e1) + ".dat").c_str(), "w");
    for (std::size_t l = 0; l < N; l++) {
        const std::size_t j = calc_j(l, nx);
        const std::size_t i = calc_i(l, j, nx);
        fprintf(file, "%ld %ld %ld %f\n", l, i, j, b[l]);
    }
    fclose(file);
}