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
