#include <array>
#include <cmath>
#include <iostream>

template <std::size_t x, std::size_t y>
using matrix = std::array < std::array<double, y>, x>;

constexpr double delta = 0.2;
constexpr std::size_t nx = 128;
constexpr std::size_t ny = 128;
constexpr double xmax = delta * nx;
constexpr double ymax = delta * ny;
constexpr double TOL = 1e-8;

constexpr double VB1(const double y) {
    return std::sin(M_PI * y / ymax);
}

constexpr double VB2(const double x) {
    return -std::sin(2 * M_PI * x / xmax);
}

constexpr double VB3(const double y) {
    return VB1(y);
}

constexpr double VB4(const double x) {
    return - VB2(x);
}

constexpr matrix<nx + 1, ny + 1> create_matrix_v() {
    matrix<nx + 1, ny + 1> result = { { 0 } };
    for (std::size_t j = 0; j < ny + 1; j++) {
        result[0][j] = VB1(delta * static_cast<double>(j));
        result[nx][j] = VB3(delta * static_cast<double>(j));
    }
    for (std::size_t i = 0; i < nx + 1; i++) {
        result[i][ny] = VB2(delta * static_cast<double>(i));
        result[i][0] = VB4(delta * static_cast<double>(i));
    }
    return result;
}

constexpr double nextV(const matrix<nx + 1, ny + 1>& V, const std::size_t i, const std::size_t j, const std::size_t k) {
    return 0.25 * (V[i + k][j] + V[i - k][j] + V[i][j + k] + V[i][j - k]);
}

constexpr void relax(matrix<nx + 1, ny + 1>& V, const std::size_t k){
    for(std::size_t i = k; i <= nx - 1; i += k){
        for(std::size_t j = k; j <= ny - 1; j += k){
            V[i][j] = nextV(V, i, j, k); 
        }
    }
}

constexpr double calcS(const matrix<nx + 1, ny + 1>& V, const std::size_t k) {
    double result = 0;

    for (int i = 0; i <= nx - k; i += k) {
        for (int j = 0; j <= ny - k; j += k) {
            double first_part = (V[i + k][j] - V[i][j]) / (2 * k * delta);
            double second_part = (V[i + k][j + k] - V[i][j + k]) / (2 * k * delta);
            const double first_sum = (first_part + second_part) * (first_part + second_part);
            first_part = (V[i][j + k] - V[i][j]) / (2 * k * delta);
            second_part = (V[i + k][j + k] - V[i + k][j]) / (2 * k * delta);
            const double second_sum = (first_part + second_part) * (first_part + second_part);
            result += k * k * delta * delta / 2 * (first_sum + second_sum);
        }
    }
    return result;
}

constexpr void we_want_more_density(matrix<nx + 1, ny + 1>& V, const std::size_t k) {
    for(std::size_t i = 0; i <= nx - k; i += k) { // maybe nx - k
        for(std::size_t j = 0; j <= ny - k; j += k) { // maybe ny - k
            V[i + k / 2][j + k / 2] = 0.25 * (V[i][j] + V[i + k][j] + V[i][j + k] + V[i + k][j + k]);
            if (j != ny - k) V[i + k / 2][j + k] = 0.5 * (V[i][j + k] + V[i + k][j + k]);
            if (i != nx - k) V[i + k][j + k / 2] = 0.5 * (V[i + k][j] + V[i + k][j + k]);
        }
    }
}

void write_to_file(const matrix<nx + 1, ny + 1>& V, const std::size_t k) {
    FILE* file = fopen(("V_k" + std::to_string(k) + ".dat").c_str(), "w");

    for(int i = 0; i <= nx - k; i += k) {
        for(int j = 0; j <= ny - k; j += k) {
            fprintf(file, "%f %f %f\n", i * delta, j * delta, V[i][j]);
        }
    }

    fclose(file);
}

void solution() {
    matrix<nx + 1, ny + 1> V = create_matrix_v();
    FILE* file = fopen("S.dat", "w");
    int i = 0;
    for(std::size_t k = 16; k > 0; k /= 2) {
        double Sit = calcS(V, k);
        double Sit_prev = 0.0;
        do {
            Sit_prev = Sit;
            relax(V, k);
            Sit = calcS(V, k);
            fprintf(file, "%d %f\n", i++, Sit);
        } while(std::abs((Sit - Sit_prev) / Sit_prev) > TOL);
        write_to_file(V, k);
        we_want_more_density(V, k);
    }
    fclose(file);

}

int main() {
    solution();
    return 0;
}