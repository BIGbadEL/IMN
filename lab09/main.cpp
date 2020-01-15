#include <cmath>
#include <iostream>

constexpr double delta = 1.0;
constexpr int nx = 40;
constexpr int ny = 40;
constexpr int n = (nx + 1) * (ny + 1);
constexpr double delta_t = 1;
constexpr int TA = 40;
constexpr int TB = 0;
constexpr int TC = 30;
constexpr int TD = 0;
constexpr int KB = 0.1;
constexpr int KD = 0.6;
constexpr int IT_MAX = 2000;

class gsl_matrix;
class gsl_vector;
class gsl_permutation;
enum CBLAS_TRANSPOSE_t { CblasNoTrans };
void gsl_matrix_set(gsl_matrix*, int, int, double);
void gsl_vector_set(gsl_vector*, int, int);
double gsl_vector_get(gsl_vector*, int);
gsl_matrix* gsl_matrix_calloc(int, int);
gsl_vector* gsl_vector_calloc(int);
gsl_permutation* gsl_permutation_calloc(int);
void gsl_linalg_LU_decomp(gsl_matrix*, gsl_permutation*, int*);



int gsl_blas_dgemv ( CBLAS_TRANSPOSE_t TransA , double alpha , const gsl_matrix * A , const gsl_vector * x , double beta , gsl_vector * y );

int gsl_blas_daxpy ( double alpha , const gsl_vector * x , gsl_vector * y );

int gsl_linalg_LU_solve ( const gsl_matrix * LU , const gsl_permutation * p , const gsl_vector * b , gsl_vector * x );
void gsl_matrix_free(gsl_matrix*);
void gsl_vector_free(gsl_vector*);
void gsl_permutation_free(gsl_permutation*);

constexpr int l(const int i, const int j) {
    return i + j * (nx + 1);
}

constexpr double x(const int i) {
    return i * delta;
}

constexpr double y(const int j) {
    return j * delta;
}

constexpr void fill_a(gsl_matrix* a) {
    for (int i = 1; i <= nx - 1; i++) {
        for (int j = 1; j <= ny - 1; j++) {
            gsl_matrix_set(a, l(i, j), l(i, j) - nx - 1, delta_t / (2 * std::pow(delta, 2)));
            gsl_matrix_set(a, l(i, j), l(i, j) - 1, delta_t / (2 * std::pow(delta, 2)));
            gsl_matrix_set(a, l(i, j), l(i, j) + 1, delta_t / (2 * std::pow(delta, 2)));
            gsl_matrix_set(a, l(i, j), l(i, j) + nx + 1, delta_t / (2 * std::pow(delta, 2)));
            gsl_matrix_set(a, l(i, j), l(i, j), -(2 * delta_t) / std::pow(delta, 2) - 1);
        }
    }

    for (int j = 0; j <= ny; j++) {
        gsl_matrix_set(a, l(0, j), l(0, j), 1);
    }

    for (int j = 0; j <= ny; j++) {
        gsl_matrix_set(a, l(nx, j), l(nx, j), 1);
    }

    for (int i = 1; i <= nx - 1; i++) {
        gsl_matrix_set(a, l(i, ny), l(i, ny) - nx - 1, -1 / (KB * delta));
        gsl_matrix_set(a, l(i, ny), l(i, ny), 1 + 1 / (KB * delta));
    }

    for (int i = 1; i <= nx - 1; i++) {
        gsl_matrix_set(a, l(i, 0), l(i, 0), 1 + 1 / (KD * delta));
        gsl_matrix_set(a, l(i, 0), l(i, 0) + nx + 1, -1 / (KD * delta));
    }
}

constexpr void fill_b(gsl_matrix* b) {
    for (int i = 1; i <= nx - 1; i++) {
        for (int j = 1; j <= ny - 1; j++) {
            gsl_matrix_set(b, l(i, j), l(i, j) - nx - 1, -delta_t / (2 * std::pow(delta, 2)));
            gsl_matrix_set(b, l(i, j), l(i, j) - 1, -delta_t / (2 * std::pow(delta, 2)));
            gsl_matrix_set(b, l(i, j), l(i, j) + 1, -delta_t / (2 * std::pow(delta, 2)));
            gsl_matrix_set(b, l(i, j), l(i, j) + nx + 1, -delta_t / (2 * std::pow(delta, 2)));
            gsl_matrix_set(b, l(i, j), l(i, j), (2 * delta_t) / std::pow(delta, 2) - 1);
        }
    }

    for (int j = 0; j <= ny; j++) {
        gsl_matrix_set(b, l(0, j), l(0, j), 1);
    }

    for (int j = 0; j <= ny; j++) {
        gsl_matrix_set(b, l(nx, j), l(nx, j), 1);
    }

    for (int k = 0; k < n; k++) {
        for (int i = 1; i <= nx - 1; i++) {
            gsl_matrix_set(b, l(i, ny), k, 0);
        }
    }

    for (int k = 0; k < n; k++) {
        for (int i = 1; i <= nx - 1; i++) {
            gsl_matrix_set(b, l(i, 0), k, 0);
        }
    }
}

constexpr void fill_c(gsl_vector* c) {

    for (int j = 0; j <= ny; j++) {
        gsl_vector_set(c, l(0, j), 0);
    }

    for (int j = 0; j <= ny; j++) {
        gsl_vector_set(c, l(nx, j), 0);
    }

    for (int i = 1; i <= nx - 1; i++) {
        gsl_vector_set(c, l(i, ny), TB);
    }

    for (int i = 1; i <= nx - 1; i++) {
        gsl_vector_set(c, l(i, 0), TD);
    }
}

constexpr void fill_T(gsl_vector* T) {
    for (int j = 0; j <= ny; j++) {
        gsl_vector_set(T, l(0, j), TA);
    }

    for (int j = 0; j <= ny; j++) {
        gsl_vector_set(T, l(nx, j), TC);
    }

    for (int i = 1; i <= nx - 1; i++) {
        for (int j = 0; j <= ny; j++) {
            gsl_vector_set(T, l(i, j), 0);
        }
    }
}

double nabla_kwadrat_razy_T(gsl_vector* T, int k) {
    return ((gsl_vector_get(T, k + 1) - 2 * gsl_vector_get(T, k) + gsl_vector_get(T, k - 1)) / std::pow(delta, 2))
        + ((gsl_vector_get(T, k + nx + 1) - 2 * gsl_vector_get(T, k) + gsl_vector_get(T, k - nx - 1)) / std::pow(delta, 2));
}

int main() {
    gsl_matrix *a = gsl_matrix_calloc(n, n);
    gsl_matrix *b = gsl_matrix_calloc(n, n);
    gsl_vector *c = gsl_vector_calloc(n);
    gsl_vector *d = gsl_vector_calloc(n);
    gsl_vector *T = gsl_vector_calloc(n);
    gsl_permutation* p = gsl_permutation_calloc(n);
    int signum = 0;

    fill_a(a);
    fill_b(b);
    fill_c(c);
    fill_T(T);


    gsl_linalg_LU_decomp(a, p, &signum);

    for( int it=1;  it<=IT_MAX; it++ ){
        gsl_blas_dgemv( CblasNoTrans, 1, b, T, 0, d );
        gsl_blas_daxpy( 1, c, d );
        gsl_linalg_LU_solve( a, p, d, T );
        if( it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000 ){
            std::cout << "creating file: out_T_" + std::to_string(it) + ".dat" << std::endl;
            FILE* f1 = fopen(("out_T_" + std::to_string(it) + ".dat").c_str(), "w");
            std::cout << "creating file: out_d2T_" + std::to_string(it) + ".dat" << std::endl;
            FILE* f2 = fopen(("out_d2T_" + std::to_string(it) + ".dat").c_str(), "w");
            for( int i=1; i<=nx-1; i++ ){
                for( int j=1; j<=ny-1; j++ ){
                    fprintf(f1, "%f %f %f\n", x(i), y(j), gsl_vector_get(T, l(i, j)));
                    fprintf(f2, "%f %f %f\n", x(i), y(j), nabla_kwadrat_razy_T(T, l(i, j)));
                }
            }
            fclose(f1);
            fclose(f2);
        }
    }

    gsl_matrix_free(a);
    gsl_matrix_free(b);
    gsl_vector_free(c);
    gsl_vector_free(T);
    gsl_permutation_free(p);
    return 0;
}