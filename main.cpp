#include <iostream>
#include <vector>
#include <cmath>

struct Point {
    double x;
    double y;
};

std::ostream& operator<<(std::ostream& stream, const std::vector<Point>& vec) {
    for(auto elem : vec){
        stream << elem.x << " " << elem.y << "\n";
    }
    return stream;
}

void vec_to_file(const std::vector<Point>& vec, const char* name){
    FILE* file = fopen(name, "wb");
    for(auto elem : vec){
        fprintf(file, "%f %f\n", elem.x, elem.y);
    }
    fclose(file);
}

constexpr double y_0 = 1;
constexpr double Q0 = 0;
constexpr double I0 = 0;
constexpr double R = 100;
constexpr double L = 0.1;
constexpr double C = 0.001;



std::vector<Point> first_euler_method(double dt, double lambda, double start, double stop) {
    std::vector<Point> result;
    result.push_back(Point{start, y_0});
    for(double i = start + dt; i <= stop; i += dt){
        Point& temp = result.back();
        result.push_back(Point{i, temp.y + dt * lambda * temp.y});
    }
    return result;
}

std::vector<Point> RK2(double dt, double lambda, double start, double stop) {
    std::vector<Point> result;
    result.push_back(Point{start, y_0});
    for(double i = start + dt; i <= stop; i += dt){
        Point& temp = result.back();
        double k1 = lambda * temp.y;
        double k2 = lambda * (temp.y + dt * k1);
        result.push_back(Point{i,temp.y + dt / 2  * (k1 + k2)});
    }
    return result;
}

void RK4_1_1(const double dt, const double lambda, const double start, const double wsp, std::vector<Point>& I, std::vector<Point>& Q) {
    double w0 = 1.0 / std::sqrt(L * C);
    double T0 = 2 * M_PI / w0;
    double wv = lambda * w0;
    double stop = T0 * wsp;
    auto V_t = [wv](double t) -> double {
        return 10.0 * std::sin(wv * t);
    };
    I.push_back(Point{start, I0});
    Q.push_back(Point{start, Q0});
    for(double i = start + dt; i <= stop; i += dt){
        double qn_1 = Q.back().y;
        double in_1 = I.back().y;
        double k1q = in_1;
        double k1i = V_t(i) / L - 1 / (L * C) * qn_1 - R / L * in_1;
        double k2q = in_1 + dt / 2 * k1i;
        double k2i = V_t(i + 0.5 * dt) / L - 1 / (L * C) * (qn_1 + dt / 2 * k1q) - R / L * (in_1 + dt / 2 * k1i);
        double k3q = in_1 + dt / 2 * k2i;
        double k3i = V_t(i + 0.5 * dt) / L - 1 / (L * C) * (qn_1 + dt / 2 * k2q) - R / L * (in_1 + dt / 2 * k2i);
        double k4q = in_1 + dt * k3i;
        double k4i = V_t(i + dt) - 1 / (L * C) * (qn_1 + dt * k3q) - R / L * (in_1 + dt * k3i);
        Q.push_back(Point{i,qn_1 + dt / 6  * (k1q + 2 * k2q + 2 * k3q + k4q)});
        I.push_back(Point{i, in_1 + dt / 6  * (k1i + 2 * k2i + 2 * k3i + k4i)});
    }
}

std::vector<Point> RK4(double dt, double lambda, double start, double stop) {
    std::vector<Point> result;
    result.push_back(Point{start, y_0});
    for(double i = start + dt; i <= stop; i += dt){
        Point& temp = result.back();
        double k1 = lambda * temp.y;
        double k2 = lambda * (temp.y + dt / 2 * k1);
        double k3 = lambda * (temp.y + dt / 2 * k2);
        double k4 = lambda * (temp.y + dt * k3);
        result.push_back(Point{i,temp.y + dt / 6  * (k1 + 2 * k2 + 2 * k3 + k4)});
    }
    return result;
}

std::vector<Point> y_t(double dt, double start, double stop, double lambda){
    std::vector<Point> result;
    for(double i = start; i <= stop; i += dt) {
        result.push_back(Point{i, std::exp(lambda * i)});
    }
    return result;
}

std::vector<Point> diff(std::vector<Point>& first, std::vector<Point>& second) {
    std::vector<Point> result;
    for(unsigned i = 0; i < first.size(); i++) {
        result.push_back(Point{first[i].x, first[i].y - second[i].y});
    }
    return result;
}

auto main() -> int {
    auto vec1 = first_euler_method(0.1, -1, 0, 5);
    vec_to_file(vec1, "1a.dat");
    auto vec1_exact = y_t(0.1, 0, 5, -1);
    vec_to_file(diff(vec1, vec1_exact), "1a_error.dat");
    vec1.clear();
    vec1_exact.clear();

    vec1 = first_euler_method(0.01, -1, 0, 5);
    vec_to_file(vec1, "1b.dat");
    vec1_exact = y_t(0.01, 0, 5, -1);
    vec_to_file(diff(vec1, vec1_exact), "1b_error.dat");
    vec1.clear();
    vec1_exact.clear();

    vec1 = first_euler_method(1, -1, 0, 5);
    vec_to_file(vec1, "1c.dat");
    vec1_exact = y_t(1, 0, 5, -1);
    vec_to_file(diff(vec1, vec1_exact), "1c_error.dat");
    vec1.clear();
    vec1_exact.clear();

    //
    vec1 = RK2(0.1, -1, 0, 5);
    vec_to_file(vec1, "2a.dat");
    vec1_exact = y_t(0.1, 0, 5, -1);
    vec_to_file(diff(vec1, vec1_exact), "2a_error.dat");
    vec1.clear();
    vec1_exact.clear();
    //

    //
    vec1 = RK2(0.01, -1, 0, 5);
    vec_to_file(vec1, "2b.dat");
    vec1_exact = y_t(0.01, 0, 5, -1);
    vec_to_file(diff(vec1, vec1_exact), "2b_error.dat");
    vec1.clear();
    vec1_exact.clear();
    //

    //vec_to_file(RK2(0.01, -1, 0, 5), "2b.dat");

    //
    vec1 = RK2(1, -1, 0, 5);
    vec_to_file(vec1, "2c.dat");
    vec1_exact = y_t(1, 0, 5, -1);
    vec_to_file(diff(vec1, vec1_exact), "2c_error.dat");
    vec1.clear();
    vec1_exact.clear();
    //
    
    //
    vec1 = RK4(0.1, -1, 0, 5);
    vec_to_file(vec1, "3a.dat");
    vec1_exact = y_t(0.1, 0, 5, -1);
    vec_to_file(diff(vec1, vec1_exact), "3a_error.dat");
    vec1.clear();
    vec1_exact.clear();
    //
    //vec_to_file(RK4(0.1, -1, 0, 5), "3a.dat");
    //

    //
    vec1 = RK4(0.01, -1, 0, 5);
    vec_to_file(vec1, "3b.dat");
    vec1_exact = y_t(0.01, 0, 5, -1);
    vec_to_file(diff(vec1, vec1_exact), "3b_error.dat");
    vec1.clear();
    vec1_exact.clear();
    //
    //vec_to_file(RK4(0.01, -1, 0, 5), "3b.dat");
    //

    //
    vec1 = RK4(1, -1, 0, 5);
    vec_to_file(vec1, "3c.dat");
    vec1_exact = y_t(1, 0, 5, -1);
    vec_to_file(diff(vec1, vec1_exact), "3c_error.dat");
    vec1.clear();
    vec1_exact.clear();
    //
    //vec_to_file(RK4(1, -1, 0, 5), "3c.dat");
    //
    
    
    std::vector<Point> Q;
    std::vector<Point> I;
    RK4_1_1(1e-4, 0.5, 0, 4, I, Q);
    vec_to_file(Q, "4aq.dat");
    vec_to_file(I, "4ai.dat");
    Q.clear();
    I.clear();
    RK4_1_1(1e-4, 0.8, 0, 4, I, Q);
    vec_to_file(Q, "4bq.dat");
    vec_to_file(I, "4bi.dat");
    Q.clear();
    I.clear();
    RK4_1_1(1e-4, 1.0, 0, 4, I, Q);
    vec_to_file(Q, "4cq.dat");
    vec_to_file(I, "4ci.dat");
    Q.clear();
    I.clear();
    RK4_1_1(1e-4, 1.2, 0, 4, I, Q);
    vec_to_file(Q, "4dq.dat");
    vec_to_file(I, "4di.dat");
    Q.clear();
    I.clear();

}