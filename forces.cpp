#include "forces.h"
#include <math.h>

//REVOIR
double LJ_potential(particle p_1, particle p_2, double L){
  double r1_2_6_i = pow(1/r(p_1, p_2, L) , 3);
  return 4*(pow(r1_2_6_i, 2) - r1_2_6_i);
}


double d_LJ_potential(particle p_1, particle p_2, double L){
  double r1_2 = r(p_1, p_2, L);
  double r1_2_i = 1./r1_2;
  double r1_2_6_i = pow(r1_2_i,3);
  double unnormed_potential = 48. * r1_2_i * r1_2_6_i * (r1_2_6_i-0.5);
  return unnormed_potential;
}


std::vector<std::vector<double>> dV_update(std::vector<particle> Particles, double L, double rc=1.5){
  std::vector<std::vector<double>> dV(Particles.size());
  for (int i=0; i<Particles.size(); i++){
    dV[i].assign(3, 0.);
  }
  double r = 0;
  for (int i = 0; i<Particles.size()-1; i++){
    for (int j = i+1; j<Particles.size(); j++){
      double d_LJ = d_LJ_potential(Particles[i], Particles[j], L);
      std::vector<double> U = direction(Particles[i], Particles[j], L);
      r = L2_norm(U);
      if (r <= rc){
        dV[i][0] -= d_LJ * U[0];
        dV[i][1] -= d_LJ * U[1];
        dV[i][2] -= d_LJ * U[2];
        dV[j][0] += d_LJ * U[0];
        dV[j][1] += d_LJ * U[1];
        dV[j][2] += d_LJ * U[2];
      }

    }
  }
  return dV;
}


void Stormer_Verlet(std::vector<particle>& Particles, double dt, double L){
  std::vector<std::vector<double>> dV = dV_update(Particles, L);
  for (int i = 0; i<Particles.size(); i++){
    Particles[i].px -= dt/2. * dV[i][0];
    Particles[i].py -= dt/2. * dV[i][1];
    Particles[i].pz -= dt/2. * dV[i][2];
    Particles[i].evolve(L, dt);
  }
  dV = dV_update(Particles, L);
  for (int i = 0; i<Particles.size(); i++){
    Particles[i].px -= dt/2. * dV[i][0];
    Particles[i].py -= dt/2. * dV[i][1];
    Particles[i].pz -= dt/2. * dV[i][2];
  }
}


void FD_term(std::vector<particle>& Particles, double dt, double gamma, double beta){
  double alpha_1, alpha_2;
  std::normal_distribution<double> G(0.0,1.0);
  std::default_random_engine generator;
  for (int i = 0; i<Particles.size(); i++){
    alpha_1 = exp(-gamma * dt);
    alpha_2 = pow(alpha_1, 2);
    Particles[i].px = alpha_1 * Particles[i].px + pow((1-alpha_2)/beta,0.5)*G(generator);
    Particles[i].py = alpha_1 * Particles[i].py + pow((1-alpha_2)/beta, 0.5)*G(generator);
    Particles[i].pz = alpha_1 * Particles[i].pz + pow((1-alpha_2)/beta, 0.5)*G(generator);
  }
}


double V(std::vector<particle> Particles, double L){
 double V = 0.;
 for (int i = 0; i < Particles.size(); i++){
   particle particle_i = Particles[i];
   for (int j = i+1; j < Particles.size() ; j++){
     particle particle_j = Particles[j];
     V += LJ_potential(particle_i, particle_j, L);
   }
 }
 return V;
}


double E(std::vector<particle> Particles, double L){
  double E = 0.;
  particle p;
  for (int i = 0; i < Particles.size(); i++){
    p = Particles[i];
    E += (pow(p.px,2)+pow(p.py,2)+pow(p.pz,2)) / (2*p.mass);
  }
  return E+V(Particles, L);
}

double T_cin(std::vector<particle> Particles, double kb){
  double u2_mean = 0;
  for (int i = 0; i<Particles.size(); i++){
    u2_mean += (pow(Particles[i].px, 2) + pow(Particles[i].py, 2) + pow(Particles[i].pz, 2)) / (Particles[i].mass);
  }
  u2_mean = u2_mean/Particles.size();
  double T = u2_mean/3;
  return T;
}
