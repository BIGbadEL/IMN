#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

constexpr double Beta = 0.001;
constexpr int N = 500;
constexpr double my_gamma = 0.1;
constexpr int tmax = 100;
constexpr double dt = 0.1;
constexpr double u0 = 1;
constexpr double Tol = 1e-6;
constexpr double alfa = Beta * N - my_gamma;
constexpr double c1 = 0.5 - 1.732050808 / 6;
constexpr double c2 = 0.5 + 1.732050808 / 6;
constexpr double a11 = 0.25;
constexpr double a12 = 0.25 - 1.732050808 / 6;
constexpr double a21 = 0.25 + 1.732050808 / 6;
constexpr double a22 = 0.25;
constexpr double b = 0.5;

struct Point {
  double x;
  double y;
};

void vec_to_file(const std::vector<Point> &vec, const char *name) {
  FILE *file = fopen(name, "wb");
  for (auto elem : vec) {
    fprintf(file, "%f %f\n", elem.x, elem.y);
  }
  fclose(file);
}

void solution(std::vector<Point> &u, std::vector<Point> &z,
              std::function<double(double, double)> fun) {
  u.push_back(Point{0.0, u0});
  for (double t = dt; t <= tmax; t += dt) {
    double u_next = u.back().y;
    int count = 0;
    double temp;
    do {
      temp = u_next;
      u_next = fun(u.back().y, u_next);
      count++;
    } while (std::abs(u_next - temp) > Tol && count <= 20);
    u.push_back(Point{t, u_next});
  }
  for (auto el : u) {
    z.push_back(Point{el.x, N - el.y});
  }
}

void RK2(std::vector<Point> &u, std::vector<Point> &z) {
  u.push_back(Point{0.0, u0});
  auto F1 = [](double u1, double u2, double un) {
    return u1 - un -
           dt * (a11 * (alfa * u1 - Beta * u1 * u1) +
                 a12 * (alfa * u2 - Beta * u2 * u2));
  };
  auto F2 = [&F1](double u2, double u1, double un) {
    return F1(u2, u1, un);
  };
  auto f = [](double t, double u) {
    return (Beta * N - my_gamma) * u - Beta * u * u;
  };

  for (double t = dt; t <= tmax; t += dt) {
    double u_next = u.back().y;
    double U1 = u_next;
    double U2 = u_next;
    double m11 = 1 - dt * a11 * (alfa - 2 * Beta * U1);
    double m12 = - dt * a12 * (alfa - 2 * Beta * U2);
    double m21 = - dt * a21 * (alfa - 2 * Beta * U1);
    double m22 = 1 - dt * a22 * (alfa - 2 * Beta * U2);

    int count = 0;
    double dU1;
    double dU2;
    do {
      dU1 = (F2(U1, U2, u_next) * m12 - F1(U1, U2, u_next) * m22) / (m11 * m22 - m12 * m21);
      dU2 = (F1(U1, U2, u_next) - F2(U1, U2, u_next)) / (m11 * m22 - m12 * m21);
      U1 += dU1;
      U2 += dU2;
      count++;
    } while (std::abs(dU1) > Tol && std::abs(dU2) > Tol && count <= 20);
    u.push_back(Point{t, u_next + dt * (b * f(t + c1 * dt, U1) + b * f(t + c2 * dt, U2))});
  }
  for (auto el : u) {
    z.push_back(Point{el.x, N - el.y});
  }
}

int main() {
  std::vector<Point> u;
  std::vector<Point> z;
  auto picard_lambda = [](double un, double u_mid) {
    return un + dt / 2 *
                    ((alfa * un - Beta * un * un) +
                     (alfa * u_mid - Beta * u_mid * u_mid));
  };
  solution(u, z, picard_lambda);
  vec_to_file(u, "u(t).dat");
  vec_to_file(z, "z(t).dat");
  u.clear();
  z.clear();
  auto newton_lambda = [](double un, double u_mid) {
    return u_mid - (u_mid - un -
                    dt / 2. *
                        ((alfa * un - Beta * un * un) +
                         (alfa * u_mid - Beta * u_mid * u_mid))) /
                       (1 - dt / 2.0 * (alfa - 2 * Beta * u_mid));
  };
  solution(u, z, newton_lambda);
  vec_to_file(u, "u(t)_new.dat");
  vec_to_file(z, "z(t)_new.dat");
  u.clear();
  z.clear();
  RK2(u, z);
  vec_to_file(u, "u(t)_RK2.dat");
  vec_to_file(z, "z(t)_RK2.dat");
  return 0;
}