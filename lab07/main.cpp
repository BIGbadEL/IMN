#include <array>
#include <cmath>
#include <iostream>
#include <string>
template <std::size_t x, std::size_t y>
using matrix = std::array<std::array<double, y>, x>;


constexpr double delta = 0.01;
constexpr double p = 1.0;
constexpr double mi = 1.0;
constexpr std::size_t nx = 200;
constexpr std::size_t ny = 90;
constexpr std::size_t i1 = 50;
constexpr std::size_t my_j1 = 55;
constexpr std::size_t j2 = my_j1 + 2;
constexpr std::size_t it_max = 20000;

constexpr double y(const std::size_t i) {
    return delta * i;
}

constexpr double x(const std::size_t i) {
    return delta * i;
}

constexpr void psi_A(matrix<nx + 1, ny + 1>& psi, const double Qwe) {
    for (std::size_t j = my_j1; j < ny + 1; j++) {
        psi[0][j] = Qwe / (2 * mi) * (std::pow(y(j), 3) / 3 - std::pow(y(j), 2) / 2 * (y(my_j1) + y(ny)) + y(j) * y(my_j1) * y(ny));
    }
}

constexpr void psi_C(matrix<nx + 1, ny + 1>& psi, const double Qwy, const double Qwe) {
    for (std::size_t j = 0; j < ny + 1; j++) {
        //psi[nx][j] = Qwe / (2 * mi) * (std::pow(y(j), 3) / 3 - std::pow(y(j), 2) / 2 * y(ny)) + (Qwe * y(my_j1) * y(my_j1) * (-y(my_j1) + 3 * y(ny))) / (12 * mi);
        psi[nx][j] = Qwy /(2*mi) * (y(j)*y(j)*y(j)/3 - y(j)*y(j)/2 * y(ny)) + Qwe *y(my_j1)*y(my_j1) * (-y(my_j1) + 3*y(ny)) / (12*mi)  ;
    }
}

constexpr void psi_B(matrix<nx + 1, ny + 1>& psi) {
    for (std::size_t i = 1; i < nx; i++) {
        psi[i][ny] = psi[0][ny];
    }
}

constexpr void psi_D(matrix<nx + 1, ny + 1>& psi) {
    for (std::size_t i = i1; i < nx; i++) {
        psi[i][0] = psi[0][my_j1];
    }
}

constexpr void psi_E(matrix<nx + 1, ny + 1>& psi) {
    for (std::size_t j = 1; j < my_j1 + 1; j++) {
        psi[i1][j] = psi[0][my_j1];
    }
}

constexpr void psi_F(matrix<nx + 1, ny + 1>& psi) {
    for (std::size_t i = 1; i < i1 + 1; i++) {
        psi[i][my_j1] = psi[0][my_j1];
    }
}

constexpr matrix<nx + 1, ny + 1> create_psi(const double Qwe, const double Qwy) {
    matrix<nx + 1, ny + 1> psi = { { 0 } };
    psi_A(psi, Qwe);
    psi_C(psi, Qwy, Qwe);
    psi_B(psi);
    psi_D(psi);
    psi_E(psi);
    psi_F(psi);
    return psi;
}

constexpr double new_psi_i_j(const std::size_t i, const std::size_t j, const matrix<nx + 1, ny + 1>& psi, const matrix<nx + 1, ny + 1>& dzeta) {
    return 0.25 * (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - delta * delta * dzeta[i][j]);
}

constexpr void dzeta_A(matrix<nx + 1, ny + 1>& dzeta, const double Qwe) {
    for (std::size_t j = my_j1; j < ny + 1; j++) {
        dzeta[0][j] = Qwe / (2 * mi) * (2 * y(j) - y(my_j1) - y(ny));
    }
}

constexpr void dzeta_C(matrix<nx + 1, ny + 1>& dzeta, const double Qwe) {
    for (std::size_t j = 0; j < ny + 1; j++) {
        dzeta[nx][j] = Qwe / (2 * mi) * (2 * y(j) - y(ny));
    }
}

constexpr void dzeta_B(matrix<nx + 1, ny + 1>& dzeta, const matrix<nx + 1, ny + 1>& psi) {
    for (std::size_t i = 1; i < nx; i++) {
        dzeta[i][ny] = 2.0 / (delta * delta) * (psi[i][ny - 1] - psi[i][ny]);
    }
}

constexpr void dzeta_D(matrix<nx + 1, ny + 1>& dzeta, const matrix<nx + 1, ny + 1>& psi) {
    for (std::size_t i = i1 + 1; i < nx; i++) {
        dzeta[i][0] = 2.0 / (delta * delta) * (psi[i][1] - psi[i][0]);
    }
}

constexpr void dzeta_E(matrix<nx + 1, ny + 1>& dzeta, const matrix<nx + 1, ny + 1>& psi) {
    for (std::size_t j = 1; j < my_j1; j++) {
        dzeta[i1][j] = 2.0 / (delta * delta) * (psi[i1 + 1][j] - psi[i1][j]);
    }
}

constexpr void dzeta_F(matrix<nx + 1, ny + 1>& dzeta, const matrix<nx + 1, ny + 1>& psi) {
    for (std::size_t i = 1; i < i1 + 1; i++) {
        dzeta[i][my_j1] = 2.0 / (delta * delta) * (psi[i][my_j1 + 1] - psi[i][my_j1]);
    }

    dzeta[i1][my_j1] = 0.5 * (dzeta[i1 - 1][my_j1] + dzeta[i1][my_j1 - 1]);
}

constexpr matrix<nx + 1, ny + 1> create_dzeta(const double Qwe, const double Qwy, const matrix<nx + 1, ny + 1>& psi) {
    matrix<nx + 1, ny + 1> dzeta = { { 0 } };
    dzeta_A(dzeta, Qwe);
    dzeta_C(dzeta, Qwy);
    dzeta_B(dzeta, psi);
    dzeta_D(dzeta, psi);
    dzeta_E(dzeta, psi);
    dzeta_F(dzeta, psi);
    return dzeta;
}

constexpr void update_dzeta(const double Qwe, const double Qwy, const matrix<nx + 1, ny + 1>& psi, matrix<nx + 1, ny + 1>& dzeta) {
    dzeta_A(dzeta, Qwe);
    dzeta_C(dzeta, Qwy);
    dzeta_B(dzeta, psi);
    dzeta_D(dzeta, psi);
    dzeta_E(dzeta, psi);
    dzeta_F(dzeta, psi);
}

constexpr double error(const matrix<nx + 1, ny + 1>& dzeta, const matrix<nx + 1, ny + 1>& psi) {
    double result = 0.0;
    for (std::size_t i = 1; i < nx; i++) {
        result += psi[i+1][j2] + psi[i-1][j2] + psi[i][j2+1] + psi[i][j2-1] - 4*psi[i][j2] - delta*delta*dzeta[i][j2];
    }
    return result;
}

constexpr double new_dzeta_i_j(const double omega, const std::size_t i, const std::size_t j, const matrix<nx + 1, ny + 1>& dzeta, const matrix<nx + 1, ny + 1>& psi) {
    //return 0.25 * (dzeta[i + 1][j] + dzeta[i - 1][j + 1] + dzeta[i][j - 1] - delta * delta * dzeta[i][j])
    //    - omega * p / (16 * mi) * ((psi[i][j + 1] - psi[i][j - 1]) * (dzeta[i + 1][j] - dzeta[i - 1][j]) - (psi[i + 1][j] - psi[i - 1][j]) * (dzeta[i][j + 1] - dzeta[i][j - 1]));
    return 0.25 * (dzeta[i+1][j] + dzeta[i-1][j] +dzeta[i][j+1] + dzeta[i][j-1]) - omega*p/(16*mi) * ((psi[i][j+1] - psi[i][j-1])*(dzeta[i+1][j] - dzeta[i-1][j])  -  (psi[i+1][j] - psi[i-1][j])*(dzeta[i][j+1] - dzeta[i][j-1])); 
}

void solve(const double Qwe) {
    const double Qwy = Qwe * (pow(y(ny), 3) - pow(y(my_j1), 3) - 3*pow(y(ny), 2)*y(my_j1) + 3*pow(y(my_j1),2)*y(ny)) / (pow(y(ny), 3));
    matrix<nx + 1, ny + 1> psi = create_psi(Qwe, Qwy);
    matrix<nx + 1, ny + 1> dzeta = {{0}};//create_dzeta(Qwe, Qwy, psi);
    FILE* error_file = fopen("error.dat", "w");
    for (std::size_t it = 1; it <= it_max; it++) {
        double omega = 1.0;
        if (it < 2000) {
            omega = 0.0;
        } else {
            omega = 1.0;
        }
        for (std::size_t i = 1; i < nx; i++) {
            for (std::size_t j = 1; j < ny; j++) {
                if (!((i <= i1 && j <= my_j1))) {
                    psi[i][j] = new_psi_i_j(i, j, psi, dzeta);
                    dzeta[i][j] = new_dzeta_i_j(omega, i, j, dzeta, psi);
                }
            }
        }
        update_dzeta(Qwe, Qwy, psi, dzeta);
        fprintf(error_file, "%f\n", error(dzeta, psi));
    }
    fclose(error_file);
    matrix<nx + 1, ny + 1> u = {{0}};
    matrix<nx + 1, ny + 1> v = {{0}};
    for (std::size_t i = 1; i < nx; i++) {
            for (std::size_t j = 1; j < ny; j++) {
                if (!((i <= i1 && j <= my_j1))) {
                    u[i][j] = (psi[i][j+1] - psi[i][j-1])/(2*delta);
                    v[i][j] = -(psi[i+1][j] - psi[i-1][j])/(2*delta);
                }
            }
        }

    FILE* out_put = fopen(("n_s_" + std::to_string(Qwe) + ".dat").c_str(), "w");
    for(std::size_t i = 0; i < nx + 1; i++){
        for(std::size_t j = 0; j < ny + 1; j++){
            fprintf(out_put, "%f %f %f %f %f %f\n", x(i), y(j), psi[i][j], dzeta[i][j], u[i][j], v[i][j]);
        }
    }
    fclose(out_put);
}

int main() {
    solve(-1000);
    solve(-4000);
    solve(4000);
    return 0;
}