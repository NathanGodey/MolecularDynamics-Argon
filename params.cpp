#include "params.h"
#include <cmath>
#include <iostream>

params::params(){
  std::cout <<std::endl <<"----SI UNITS----" <<std::endl;
  std::cout <<"TEMPERATURE : "<<this->RAW_TEMPERATURE <<" K"<<std::endl;
  std::cout <<"DENSITY : "<<this->RAW_DENSITY <<" kg/m³"<<std::endl;
  this->NB_PARTICLES = pow(this->NB_PART_PER_DIM, 3);
  this->A = pow(this->RAW_ARGON_MASS / this->RAW_DENSITY, 1/3.);
  std::cout <<"NB_PARTICLES : "<<this->NB_PARTICLES <<std::endl;
  this->L = this->NB_PART_PER_DIM * this->A;
  std::cout <<"L : "<<this->L <<" m"<<std::endl;
  std::cout <<"A : "<<this->A <<" m" <<std::endl;
  std::cout <<"dt : "<<this->RAW_SCHEME_DT <<" s" <<std::endl;
  this->L /= pow(2, 1/6.) * this->RAW_D_ARGON;
  this->A /= pow(2, 1/6.) * this->RAW_D_ARGON;


  this->K_BOLTZ = 1./119.8;
  this->TEMPERATURE = this->RAW_TEMPERATURE / 119.8;
  this->BETA = 1. / (this->TEMPERATURE);


  this->D_ARGON = 1./pow(2, 1/6.);

  this->DURATION = this->NB_ITER * this->SCHEME_DT;
  this->ARGON_MASS = 1;

  this->INI_P_SCALE = pow(3*this->TEMPERATURE * this->ARGON_MASS, 0.5);
  this->E0 = 1.;
  this->SCHEME_DT = this->RAW_SCHEME_DT/(this->RAW_D_ARGON * pow(this->RAW_ARGON_MASS / this->RAW_E_0, 0.5));



  this->DENSITY = this->ARGON_MASS * this->NB_PARTICLES / pow(this->L, 3);


  std::cout <<std::endl <<"----REDUCED UNITS----" <<std::endl;
  std::cout <<"L : "<<this->L <<std::endl;
  std::cout <<"A : "<<this->A <<std::endl;
  std::cout <<"TEMPERATURE : "<<this->TEMPERATURE <<std::endl;
  std::cout <<"DENSITY : "<<this->DENSITY <<std::endl;
  std::cout <<"DELTA_T : "<<this->SCHEME_DT <<std::endl;
  std::cout <<"BETA : " <<this->BETA <<std::endl;
  std::cout <<"GAMMA : " <<this->GAMMA <<std::endl;

}
