#include "params.h"
#include <cmath>
#include <iostream>

params::params(){
  this->NB_PARTICLES = floor(this->RAW_DENSITY * pow(this->RAW_L, 3) / this->RAW_ARGON_MASS);
  std::cout <<"NB_PARTICLES : "<<this->NB_PARTICLES <<std::endl;
  this->L = this->RAW_L/(pow(2, 1./6.)*this->RAW_D_ARGON);
  std::cout <<"L : "<<this->L <<std::endl;
  this->K_BOLTZ = 1./119.8;
  this->TEMPERATURE = this->RAW_TEMPERATURE / 119.8;
  std::cout <<"TEMPERATURE : "<<this->TEMPERATURE <<std::endl;
  this->BETA = 1. / (this->K_BOLTZ * this->TEMPERATURE);
  this->SCHEME_DT = this->RAW_SCHEME_DT * pow(10., 12);
  this->D_ARGON = 1.;
  std::cout <<"D_ARGON : "<<this->D_ARGON <<std::endl;
  this->DURATION = this->NB_ITER * this->SCHEME_DT;
  this->ARGON_MASS = this->RAW_ARGON_MASS / this->U;
  this->A = this->L / pow(this->NB_PARTICLES+1, 1./this->DIMENSION);
  this->INI_P_SCALE = pow(this->K_BOLTZ * this->TEMPERATURE * this->ARGON_MASS);
  this->E0 = 1.;
  std::cout <<"A : "<<this->A <<std::endl;
}
