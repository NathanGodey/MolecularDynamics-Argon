#pragma once
#include <vector>
#include <random>

class particle {
public:
  double x, y, z, px, py, pz;
  double mass;
  particle();
  particle(double x, double y, double z, double px, double py, double pz);
  void evolve(double L, double dt);
};

std::vector<particle> init_Argon(double A, int NB_PART_PER_DIM, double INI_P_SCALE, double MASS);
std::vector<double> random_vector(double norm, int dimension);
double L2_norm(std::vector<double> u);
std::vector<double> direction(particle p_1, particle p_2, double L);
double r(particle p_1, particle p_2, double L);
double mean(std::vector<double> sample);
double std_dev(std::vector<double> sample);
