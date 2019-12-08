#include <iostream>
#include <vector>
#include <memory>
#include <vector>
#include "fill_vec.h"
#include "mgmres.h"

void solve(std::size_t nx, std::size_t ny, int nz_num, std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& a, std::vector<double>& b, std::vector<double>& V){
    static int count = 1;
    const std::size_t N = (nx + 1) * (ny + 1);
    V.clear();
    V.resize(N);
    int iter_max = 500;
    int mr = 500;
    double tol_abs = 1e-8;
    double tol_rel = 1e-8;
    pmgmres_ilu_cr(N, nz_num, ia.data(), ja.data(), a.data(), V.data(), b.data(), iter_max, mr, tol_abs, tol_rel);
    FILE* file = fopen(("V_" + std::to_string(count) + ".dat").c_str(), "w");
    for(std::size_t l = 0; l < N; l++){
        const std::size_t j = calc_j(l, nx);
        const std::size_t i = calc_i(l, j, nx);
        fprintf(file, "%f %f %f\n", i * delta, j * delta, V[l]);
    }
    fclose(file);
    count++;
}


int main() {
    std::vector<int> ia;
    std::vector<double> b;
    std::vector<int> ja;
    std::vector<double> a;
    std::vector<double> V;
    double V1 = 10;
    double V3 = 10;
    double V2 = -10;
    double V4 = -10;
    int nz_num = 0;
    std::size_t nx = 4;
    std::size_t ny = 4;
    double e1 = 1.0;
    double e2 = 1.0;
    fill_vecs(nx, ny, e1, e2,  b, ia, a, ja, nz_num, V1, V2, V3, V4);
    solve(nx, ny, nz_num, ia, ja, a, b, V);

    nx = 50;
    ny = 50;

    fill_vecs(nx, ny, e1, e2,  b, ia, a, ja, nz_num, V1, V2, V3, V4);
    solve(nx, ny, nz_num, ia, ja, a, b, V);

    nx = 100;
    ny = 100;

    fill_vecs(nx, ny, e1, e2,  b, ia, a, ja, nz_num, V1, V2, V3, V4);
    solve(nx, ny, nz_num, ia, ja, a, b, V);

    nx = 200;
    ny = 200;

    fill_vecs(nx, ny, e1, e2,  b, ia, a, ja, nz_num, V1, V2, V3, V4);
    solve(nx, ny, nz_num, ia, ja, a, b, V);

    V1 = 0;
    V3 = 0;
    V2 = 0;
    V4 = 0;
    nx = 100;
    ny = 100;

    fill_vecs(nx, ny, e1, e2,  b, ia, a, ja, nz_num, V1, V2, V3, V4, false);
    solve(nx, ny, nz_num, ia, ja, a, b, V);

    e1 = 1.0;
    e2 = 2.0;

    fill_vecs(nx, ny, e1, e2,  b, ia, a, ja, nz_num, V1, V2, V3, V4, false);
    solve(nx, ny, nz_num, ia, ja, a, b, V);

    e1 = 1.0;
    e2 = 10.0;

    fill_vecs(nx, ny, e1, e2,  b, ia, a, ja, nz_num, V1, V2, V3, V4, false);
    solve(nx, ny, nz_num, ia, ja, a, b, V);

    return 0;
}