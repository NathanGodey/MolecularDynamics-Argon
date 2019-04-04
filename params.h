#pragma once
#include <cmath>
class params{
public:
  double RAW_L = 10. * pow(10,-9);
  double RAW_DENSITY = pow(10,3);
  int NB_PARTICLES;
  double RAW_TEMPERATURE = 300.;
  double RAW_K_BOLTZ = 1.38064852 * pow(10,-23);
  double RAW_E_0 = K_BOLTZ * 119.8;
  double RAW_D_ARGON = 0.3405 * pow(10,-9);
  double RAW_BETA = 1. / (K_BOLTZ * TEMPERATURE);
  double RAW_U = 1.6605389 * pow(10,-27);
  double RAW_ARGON_MASS = 39.948 * U;
  double RAW_SCHEME_DT = 5 * pow(10,-15);
  int NB_ITER = 100;
  int DIMENSION = 3;
  double GAMMA = 0.01;
  double INI_P_SCALE, DURATION, A, L, DENSITY, TEMPERATURE, K_BOLTZ, E0, D_ARGON, BETA, ARGON_MASS, SCHEME_DT;
  params();
};
