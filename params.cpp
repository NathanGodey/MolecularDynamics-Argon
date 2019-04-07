#include "params.h"
#include <cmath>
#include <iostream>

params::params(){

  this->NB_PARTICLES = pow(this->NB_PART_PER_DIM, 3);
  this->A = pow(this->RAW_ARGON_MASS / this->RAW_DENSITY, 1/3.);
  std::cout <<"NB_PARTICLES : "<<this->NB_PARTICLES <<std::endl;
  this->L = this->NB_PART_PER_DIM * this->A;
  std::cout <<"L : "<<this->L <<std::endl;
  std::cout <<"A : "<<this->A <<std::endl;
  this->L /= pow(2, 1/6.) * this->RAW_D_ARGON;
  this->A /= pow(2, 1/6.) * this->RAW_D_ARGON;


  this->K_BOLTZ = 1./119.8;
  this->TEMPERATURE = this->RAW_TEMPERATURE / 119.8;
  std::cout <<"TEMPERATURE : "<<this->RAW_TEMPERATURE <<std::endl;
  this->BETA = 1. / (this->K_BOLTZ * this->TEMPERATURE);
  // this->SCHEME_DT = this->RAW_SCHEME_DT / (pow(2, 1/6.) * this->RAW_D_ARGON * pow(this->RAW_ARGON_MASS / this->RAW_E_0, 0.5));
  this->SCHEME_DT = 0.5;
  this->D_ARGON = 1./pow(2, 1/6.);

  this->DURATION = this->NB_ITER * this->SCHEME_DT;
  this->ARGON_MASS = this->RAW_ARGON_MASS / this->U;

  this->INI_P_SCALE = pow(this->K_BOLTZ * this->TEMPERATURE * this->ARGON_MASS, 0.5);
  this->E0 = 1.;

  this->DENSITY = this->RAW_DENSITY * pow(2, 0.5) * pow(this->RAW_D_ARGON,3)/this->RAW_ARGON_MASS;


  std::cout <<std::endl <<"----REDUCED UNITS----" <<std::endl;
  std::cout <<"L : "<<this->L <<std::endl;
  std::cout <<"A : "<<this->A <<std::endl;
  std::cout <<"TEMPERATURE : "<<this->DENSITY <<std::endl;
  std::cout <<"DENSITY : "<<this->DENSITY <<std::endl;
  std::cout <<"DELTA_T : "<<this->SCHEME_DT <<std::endl;

}
