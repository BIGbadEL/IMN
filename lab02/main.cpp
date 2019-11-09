#include <cmath>
#include <iostream>
#include <vector>
#include <functional>

constexpr double Beta = 0.001;
constexpr int N = 500;
constexpr double my_gamma = 0.1;
constexpr int tmax = 100;
constexpr double dt = 0.1;
constexpr double u0 = 1;
constexpr double Tol = 1e-6;
constexpr double alfa = Beta * N - my_gamma;

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

void solution(std::vector<Point> &u, std::vector<Point> &z, std::function<double(double , double)> fun) {
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
                       (1 - dt / 2.0 * ( alfa - 2 * Beta * u_mid));
  };
  solution(u , z, newton_lambda);
  vec_to_file(u, "u(t)_new.dat");
  vec_to_file(z, "z(t)_new.dat");
  return 0;
}