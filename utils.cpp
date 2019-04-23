#include "utils.h"
#include <iostream>

particle::particle(double x, double y, double z, double px, double py, double pz){
  this->x = x;
  this->y = y;
  this->z = z;
  this->px = px;
  this->py = py;
  this->pz = pz;
}
particle::particle(){
}

void torify(double& param, double L){
  if (param < 0){
    param += L;
  }
  else if (param >= L){
    param -= L;
  }
  if (param<0 || param > L){
    std::cerr <<std::endl <<"Bad parameter choice : particles moving too fast relative to dimensions" <<std::endl;
    std::exit(1);
  }
}

void particle::evolve(double L, double dt){
  this->x += (dt * this->px / this->mass);
  torify(this->x, L);
  this->y += (dt * this->py / this->mass);
  torify(this->y, L);
  this->z += (dt * this->pz / this->mass);
  torify(this->z, L);
}

double minabs(double a, double b){
  if (std::abs(a) >= std::abs(b)){
    return b;
  }
  else{
    return a;
  }
}


std::vector<double> direction(particle p_1, particle p_2, double L){
  std::vector<double> U;
  double X=0, Y=0, Z=0;
  if (p_1.x <= p_2.x){
    X = minabs(p_1.x - p_2.x, p_1.x - (p_2.x-L));
  }
  else{
    X = minabs(p_1.x - p_2.x, (p_1.x-L) - p_2.x);
  }
  if (p_1.y <= p_2.y){
    Y = minabs(p_1.y - p_2.y, p_1.y - (p_2.y-L));
  }
  else{
    Y = minabs(p_1.y - p_2.y, (p_1.y-L) - p_2.y);
  }
  if (p_1.z <= p_2.z){
    Z = minabs(p_1.z - p_2.z, p_1.z - (p_2.z-L));
  }
  else{
    Z = minabs(p_1.z - p_2.z, (p_1.z-L) - p_2.z);
  }
  U.push_back(X);
  U.push_back(Y);
  U.push_back(Z);
  return U;
}

double r(particle p_1, particle p_2, double L){
  std::vector<double> U = direction(p_1, p_2, L);
  return pow(U[0],2)+pow(U[1],2)+pow(U[2],2);
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


double L2_norm(std::vector<double> u){
  double s=0;
  for (int i=0;i<u.size();i++){
    s+=pow(u[i],2);
  }
  return sqrt(s);
}


std::vector<double> random_vector(double norm, int dimension){
  std::vector<double> d;
  for (int i=0;i<dimension;i++){
    d.push_back(((double)rand()/RAND_MAX)*2 -1);
  }
  double n = L2_norm(d);
  for (int i=0;i<dimension;i++){
    d[i] *= norm / n;
  }
  return d;
}


std::vector<particle> init_Argon(double A, int NB_PART_PER_DIM, double INI_P_SCALE, double MASS){
  std::vector<particle> Particles;
  double current_X, current_Y, current_Z, PX, PY, PZ;
  for (int i = 0; i<NB_PART_PER_DIM; i++){
    current_X = i * A + A/2;
    for (int j = 0; j<NB_PART_PER_DIM; j++){
      current_Y = j * A + A/2;
      for (int k = 0; k<NB_PART_PER_DIM; k++){
        current_Z = k * A + A/2;
        std::vector<double> u = random_vector(INI_P_SCALE, 3);
        PX = u[0];
        std::cout <<"PX = "<<PX <<std::endl;
        PY = u[1];
        std::cout <<"PY = "<<PY <<std::endl;
        PZ = u[2];
        std::cout <<"PZ = "<<PZ <<std::endl;
        particle P(current_X, current_Y, current_Z, PX, PY, PZ);
        P.mass = MASS;
        Particles.push_back(P);
      }
    }
  }
  return Particles;
}
