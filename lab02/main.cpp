#include <cmath>
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

void picard(std::vector<Point> &u, std::vector<Point> &z) {
  u.push_back(Point{0.0, u0});
  auto u_n_plus_1 = [](double un, double u_mid) {
    return un + dt / 2 *
                    ((alfa * un - Beta * un * un) +
                     (alfa * u_mid - Beta * u_mid * u_mid));
  };
  for (double t = dt; t <= tmax; t += dt) {
    double u_next = u.back().y;
    int count = 0;
    double temp;
    do {
      temp = u_next;
      u_next = u_n_plus_1(u.back().y, u_next);
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
  picard(u, z);
  vec_to_file(u, "u(t).dat");
  vec_to_file(z, "z(t).dat");
  return 0;
}