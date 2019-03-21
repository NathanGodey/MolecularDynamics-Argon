#include "utils.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>
#include <fstream>
#include <random>




double L = 10. * pow(10,-7);
int NB_PARTICLES = 200;
double TEMPERATURE = 3000.;
double k_boltz = 1.38064852 * pow(10,-23);
double E_0 = k_boltz * 119.8;
double d_argon = 0.3405 * pow(10,-9);
double beta = 1. / (k_boltz * TEMPERATURE);
double u = 1.66054 * pow(10,-27);
double SCHEME_DT = 5 * pow(10,-15);



double direction_indicator(double from, double to){
  if (from > to){
    return 1.;
  }
  else if (from < to){
    return -1.;
  }
  else{
    return 0.;
  }
}


double r(particle p_1, particle p_2){
  return sqrt(pow(p_1.x-p_2.x,2)+pow(p_1.y-p_2.y,2)+pow(p_1.z-p_2.z,2));
}


double LJ_potential(particle p_1, particle p_2){
  double r1_2 = r(p_1, p_2);
  return 4. * E_0 * (pow(d_argon/r1_2, 12) - pow(d_argon/r1_2, 6));
}


double d_LJ_potential(particle p_1, particle p_2){
  double r1_2 = r(p_1, p_2);
  double unnormed_potential = -12. * pow((d_argon/r1_2),13) + 6. * pow((d_argon/r1_2),7);
  return 4. * E_0 * unnormed_potential;
}


std::vector<double> dV(std::vector<particle> Particles, int i){
  std::vector<double> dV{0., 0., 0.};
  particle particle_i = Particles[i];
  for (int j = 0; j<NB_PARTICLES && j!=i; j++){
    particle particle_j = Particles[j];
    double d_LJ = d_LJ_potential(particle_i, particle_j);
    dV[0] += d_LJ * direction_indicator(particle_i.x, particle_j.x);
    dV[1] += d_LJ * direction_indicator(particle_i.y, particle_j.y);
    dV[2] += d_LJ * direction_indicator(particle_i.z, particle_j.z);
  }
  return dV;
}

void Stormer_Verlet(std::vector<particle>& Particles, double dt = SCHEME_DT){
  std::vector<double> dV_i;
  for (int i = 0; i<NB_PARTICLES; i++){
    dV_i = dV(Particles, i);
    Particles[i].px -= dt/2. * dV_i[0];
    Particles[i].py -= dt/2. * dV_i[1];
    Particles[i].pz -= dt/2. * dV_i[2];
    Particles[i].evolve(L, dt);
    dV_i = dV(Particles, i);
    Particles[i].px -= dt/2. * dV_i[0];
    Particles[i].py -= dt/2. * dV_i[1];
    Particles[i].pz -= dt/2. * dV_i[2];
  }
}

void FD_term(std::vector<particle>& Particles, double dt = SCHEME_DT){
  double gamma, alpha_1, alpha_2;
  std::default_random_engine generator;
  std::normal_distribution<double> G(0.0,pow(dt,2));
  for (int i = 0; i<NB_PARTICLES; i++){
    gamma = 0.001 * Particles[i].mass/dt;
    alpha_1 = exp(-gamma * dt/Particles[i].mass);
    alpha_2 = pow(alpha_1, 2);
    Particles[i].px = alpha_1 * Particles[i].px + sqrt((1-alpha_2)*Particles[i].mass/beta)*G(generator);
    Particles[i].py = alpha_1 * Particles[i].py + sqrt((1-alpha_2)*Particles[i].mass/beta)*G(generator);
    Particles[i].pz = alpha_1 * Particles[i].pz + sqrt((1-alpha_2)*Particles[i].mass/beta)*G(generator);
  }
}

double V(std::vector<particle> Particles){
 double V = 0.;
 for (int i = 0; i < NB_PARTICLES; i++){
   particle particle_i = Particles[i];
   for (int j = i+1; j < NB_PARTICLES ; j++){
     particle particle_j = Particles[j];
     V += LJ_potential(particle_i, particle_j);
   }
 }
 return V;
}


double E(std::vector<particle> Particles){
  double E = 0.;
  particle p;
  for (int i = 0; i < NB_PARTICLES; i++){
    p = Particles[i];
    E += (pow(p.px,2)+pow(p.py,2)+pow(p.pz,2)) / (2*p.mass);
  }
  return E + V(Particles);
}


std::vector<particle> init_Argon(){
  std::vector<particle> Particles;
  double current_X = (L) / pow(NB_PARTICLES+1,0.33), current_Y = (L) / pow(NB_PARTICLES+1,0.33), current_Z = (L) / pow(NB_PARTICLES+1,0.33);
  for (int i = 0; i<NB_PARTICLES; i++){
    current_X += (L) / pow(NB_PARTICLES+1,0.33);
    if (current_X > L){
      current_X = fmod(current_X ,L);
      current_Y += (L) / pow(NB_PARTICLES+1,0.33);
      if (current_Y > L){
        current_Y = fmod(current_Y ,L);
        current_Z += (L) / pow(NB_PARTICLES+1,0.33);
        if (current_Z > L){
          current_Z = fmod(current_Z ,L);
        }
      }

    }
    double v_scale = pow(10., -19);
    double P_X = 5*v_scale *(((double)rand()/RAND_MAX)*2 -1),
    P_Y = 5*v_scale*(((double)rand()/RAND_MAX)*2 -1),
    P_Z = 5*v_scale*(((double)rand()/RAND_MAX)*2 -1);
    particle P(current_X, current_Y, current_Z, P_X, P_Y, P_Z);
    P.mass = 39.948 * u;
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
  int NB_ITER = 1000;
  std::ofstream pos_output;
  std::ofstream energy_output;
  pos_output.open("./positions.txt");
  energy_output.open("./energy.txt");
  pos_output <<L <<std::endl;
  pos_output <<NB_PARTICLES <<std::endl;
  pos_output <<NB_ITER <<std::endl;
  energy_output <<NB_ITER <<std::endl;
  for (int t = 0; t<NB_ITER; t++){
    std::cout <<"Itération : " <<t+1 <<"/" <<NB_ITER <<std::flush;
    for (int i=0; i<NB_PARTICLES;i++){
      particle p = Arg_Particles[i];
      pos_output <<p.x <<" " <<p.y <<" " <<p.z <<std::endl;
    }
    double current_E = E(Arg_Particles);
    E_vect.push_back(current_E);
    energy_output  <<std::setprecision(10) <<current_E <<std::endl;
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
