#include "utils.h"
#include "params.h"
#include "forces.h"
#include "test.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>
#include <fstream>
#include <random>







int main(){
  bool TEST_MODE = false;
  if (TEST_MODE){
    Test test = Test(2);
    test.exe(LJ_potential);
  }
  else{
    bool SHOW_ITER = false;
    params my_params;
    std::vector<particle> Arg_Particles = init_Argon(my_params.A, my_params.NB_PART_PER_DIM, my_params.INI_P_SCALE, my_params.ARGON_MASS);
    std::cout <<"Gaz généré" <<std::endl;
    std::vector<double> E_vect;
    std::ofstream pos_output;
    std::ofstream energy_output;
    std::ofstream tempcin_output;
    std::ofstream p_output;
    pos_output.open("./positions.txt");
    energy_output.open("./energy.txt");
    tempcin_output.open("./tempcin.txt");
    p_output.open("./momentum.txt");
    pos_output <<my_params.RAW_E_0 <<std::endl;
    pos_output <<my_params.L <<std::endl;
    pos_output <<my_params.NB_PARTICLES <<std::endl;
    pos_output <<my_params.NB_ITER <<std::endl;
    energy_output <<my_params.NB_ITER <<std::endl;
    tempcin_output <<my_params.NB_ITER <<std::endl;
    p_output <<my_params.NB_ITER <<std::endl;
    p_output <<my_params.NB_PARTICLES <<std::endl;
    for (int t = 0; t<my_params.NB_ITER; t++){
      if (SHOW_ITER) {std::cout <<"Itération : " <<t+1 <<"/" <<my_params.NB_ITER <<std::flush;}
      for (int i=0; i<my_params.NB_PARTICLES;i++){
        particle p = Arg_Particles[i];
        p_output <<L2_norm({p.px, p.py, p.pz}) <<std::endl;
        pos_output <<p.x <<" " <<p.y <<" " <<p.z <<std::endl;
      }
      double current_E = E(Arg_Particles, my_params.L);
      E_vect.push_back(current_E);
      energy_output  <<current_E <<std::endl;
      tempcin_output  <<T_cin(Arg_Particles, my_params.K_BOLTZ) <<std::endl;
      Stormer_Verlet(Arg_Particles, my_params.SCHEME_DT, my_params.L);
      FD_term(Arg_Particles, my_params.SCHEME_DT, my_params.GAMMA, my_params.BETA);
      if (SHOW_ITER) {std::cout <<"\r";}
    }
    pos_output.close();
    std::cout <<"Moyenne de E : " <<mean(E_vect) <<std::endl;
    std::cout <<"Ecart-type de E : " <<std_dev(E_vect) <<std::endl;
    energy_output <<mean(E_vect) <<std::endl;
    energy_output <<std_dev(E_vect) <<std::endl;
    energy_output.close();
    system("python3 visualizer.py");
  }

  return 0;
}
