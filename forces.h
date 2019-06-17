#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include "utils.h"


double LJ_potential(particle p_1, particle p_2, double L);
double d_LJ_potential(particle p_1, particle p_2, double L);
std::vector<std::vector<double>> dV_update(std::vector<particle> Particles);
void Stormer_Verlet(std::vector<particle>& Particles, double dt, double L);
void FD_term(std::vector<particle>& Particles, double dt, double gamma, double beta);
double E(std::vector<particle> Particles, double L);
double V(std::vector<particle> Particles, double L);
double T_cin(std::vector<particle> Particles, double kb);
