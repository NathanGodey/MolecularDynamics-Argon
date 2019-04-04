#include "utils.h"
#include "params.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>
#include <fstream>
#include <random>




params my_params;





double r(particle p_1, particle p_2){
  return sqrt(pow(p_1.x-p_2.x,2)+pow(p_1.y-p_2.y,2)+pow(p_1.z-p_2.z,2));
}


double LJ_potential(particle p_1, particle p_2){
  double r1_2 = r(p_1, p_2);
  return (pow(1./r1_2, 12) - pow(1./r1_2, 6));
}


double d_LJ_potential(particle p_1, particle p_2){
  double r1_2 = r(p_1, p_2);
  double r1_2_i = 1./r1_2;
  double r1_2_6_i = pow(r1_2_i,3);
  double unnormed_potential = 48. * r1_2_i * r1_2_6_i * (r1_2_6_i-0.5);
  return unnormed_potential;
}

std::vector<double> direction(particle p_1, particle p_2){
  std::vector<double> U;
  double X=0, Y=0, Z=0;
  if (p_1.x != p_2.x){
    X = fmod(p_1.x - p_2.x, my_params.L);
  }
  if (p_1.y != p_2.y){
    Y = fmod(p_1.x - p_2.x, my_params.L);
  }
  if (p_1.z != p_2.z){
    Z = fmod(p_1.x - p_2.x, my_params.L);
  }
  U.push_back(X);
  U.push_back(Y);
  U.push_back(Z);
  return U;
}

std::vector<std::vector<double>> dV_update(std::vector<particle> Particles){
  std::vector<std::vector<double>> dV(my_params.NB_PARTICLES);
  for (int i=0; i<my_params.NB_PARTICLES; i++){
    dV[i].assign(3, 0.);
  }
  for (int i = 0; i<my_params.NB_PARTICLES-1; i++){
    for (int j = i+1; j<my_params.NB_PARTICLES; j++){
      double d_LJ = d_LJ_potential(Particles[i], Particles[j]);
      std::vector<double> U = direction(Particles[i], Particles[j]);
      dV[i][0] += d_LJ * U[0];
      dV[i][1] += d_LJ * U[1];
      dV[i][2] += d_LJ * U[2];
      dV[j][0] -= d_LJ * U[0];
      dV[j][1] -= d_LJ * U[1];
      dV[j][2] -= d_LJ * U[2];
    }
  }
  return dV;
}

void Stormer_Verlet(std::vector<particle>& Particles, double dt = my_params.SCHEME_DT){
  std::vector<std::vector<double>> dV = dV_update(Particles);
  for (int i = 0; i<my_params.NB_PARTICLES; i++){
    Particles[i].px -= dt/2. * dV[i][0];
    Particles[i].py -= dt/2. * dV[i][1];
    Particles[i].pz -= dt/2. * dV[i][2];
    Particles[i].evolve(my_params.L, dt);
  }
  dV = dV_update(Particles);
  for (int i = 0; i<my_params.NB_PARTICLES; i++){
    Particles[i].px -= dt/2. * dV[i][0];
    Particles[i].py -= dt/2. * dV[i][1];
    Particles[i].pz -= dt/2. * dV[i][2];
  }
}

void FD_term(std::vector<particle>& Particles, double dt = my_params.SCHEME_DT){
  double alpha_1, alpha_2;
  std::default_random_engine generator;
  std::normal_distribution<double> G(0.0,1.0);
  for (int i = 0; i<my_params.NB_PARTICLES; i++){
    alpha_1 = exp(-my_params.GAMMA * dt/Particles[i].mass);
    alpha_2 = pow(alpha_1, 2);
    Particles[i].px = alpha_1 * Particles[i].px + sqrt((1-alpha_2)*Particles[i].mass/my_params.BETA)*G(generator);
    Particles[i].py = alpha_1 * Particles[i].py + sqrt((1-alpha_2)*Particles[i].mass/my_params.BETA)*G(generator);
    Particles[i].pz = alpha_1 * Particles[i].pz + sqrt((1-alpha_2)*Particles[i].mass/my_params.BETA)*G(generator);
  }
}

double V(std::vector<particle> Particles){
 double V = 0.;
 for (int i = 0; i < my_params.NB_PARTICLES; i++){
   particle particle_i = Particles[i];
   for (int j = i+1; j < my_params.NB_PARTICLES ; j++){
     particle particle_j = Particles[j];
     V += LJ_potential(particle_i, particle_j);
   }
 }
 return V;
}


double E(std::vector<particle> Particles){
  double E = 0.;
  particle p;
  for (int i = 0; i < my_params.NB_PARTICLES; i++){
    p = Particles[i];
    E += (pow(p.px,2)+pow(p.py,2)+pow(p.pz,2)) / (2*p.mass);
  }
  std::cout <<"E=" <<E <<std::endl;
  std::cout <<"V=" <<V(Particles) <<std::endl;
  return E + V(Particles);
}


std::vector<particle> init_Argon(){
  std::vector<particle> Particles;
  double current_X = my_params.A/2, current_Y = my_params.A/2, current_Z = my_params.A/2;
  for (int i = 0; i<my_params.NB_PARTICLES; i++){
    current_X += my_params.A;
    if (current_X > my_params.L){
      current_X = my_params.A/2;
      current_Y += my_params.A;
      if (current_Y > my_params.L){
        current_Y = my_params.A/2;
        current_Z += my_params.A;
        if (current_Z > my_params.L){
          current_Z = my_params.A/2;
        }
      }

    }
    double p_scale = my_params.INI_P_SCALE;
    double P_X = 5*p_scale *(((double)rand()/RAND_MAX)*2 -1),
    P_Y = 5*p_scale*(((double)rand()/RAND_MAX)*2 -1),
    P_Z = 5*p_scale*(((double)rand()/RAND_MAX)*2 -1);
    particle P(current_X, current_Y, current_Z, P_X, P_Y, P_Z);
    P.mass = my_params.ARGON_MASS;
    Particles.push_back(P);
  }
  return Particles;
}

double mean(std::vector<double> sample){
  double S =0.;
  for (int i=0; i<sample.size(); i++){
    S+=sample[i];
  }
  return S/sample.size();
}

double std_dev(std::vector<double> sample){
  double S = 0.;
  double this_mean = mean(sample);
  for (int i=0; i<sample.size(); i++){
    S+=pow(sample[i]-this_mean,2);
  }
  return sqrt(S);
}








int main(){
  std::vector<particle> Arg_Particles = init_Argon();
  std::cout <<"Gaz généré" <<std::endl;
  std::vector<double> E_vect;
  std::ofstream pos_output;
  std::ofstream energy_output;
  pos_output.open("./positions.txt");
  energy_output.open("./energy.txt");
  pos_output <<my_params.L <<std::endl;
  pos_output <<my_params.NB_PARTICLES <<std::endl;
  pos_output <<my_params.NB_ITER <<std::endl;
  energy_output <<my_params.NB_ITER <<std::endl;
  for (int t = 0; t<my_params.NB_ITER; t++){
    std::cout <<"Itération : " <<t+1 <<"/" <<my_params.NB_ITER <<std::flush;
    for (int i=0; i<my_params.NB_PARTICLES;i++){
      particle p = Arg_Particles[i];
      pos_output <<p.x <<" " <<p.y <<" " <<p.z <<std::endl;
    }
    double current_E = E(Arg_Particles);
    E_vect.push_back(current_E);
    energy_output  <<current_E <<std::endl;
    Stormer_Verlet(Arg_Particles);
    FD_term(Arg_Particles);
    std::cout <<"\r";
  }
  pos_output.close();
  std::cout <<"Moyenne de E : " <<mean(E_vect) <<std::endl;
  std::cout <<"Ecart-type de E : " <<std_dev(E_vect) <<std::endl;
  energy_output <<mean(E_vect) <<std::endl;
  energy_output <<std_dev(E_vect) <<std::endl;
  energy_output.close();
  system("python3 visualizer.py");
  return 0;
}
