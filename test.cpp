#include "params.h"
#include "forces.h"
#include "utils.h"
#include "test.h"
#include <fstream>
#define getName(var) #var

Test::Test(){
}

Test::Test(int t){
    this->type = t;
    if(t==1){std::cout << "test de superposition" << std::endl;}
    else if(t==2){std::cout << "test force" << std::endl;}
    else if(t==3){std::cout << "test potentiel" << std::endl;}
    else if(t==4){std::cout << "test rebond" << std::endl;}
    else {std::cerr << "unvalid test" << std::endl;}
}

void Test::show_rep(){
    if (this->rep==0) {std::cout << "no problem with test " << this->type << std::endl;}
    else {std::cerr << "problem with test" << this->type << std::endl;}

}

void Test::exe(std::vector<particle> particles){
    if (this->type == 1){
        int N=particles.size();
        for(int i=0;(i<N-1) && (this->rep==0);i++){
            for(int j=i+1;j<N && (this->rep==0);j++){
                if (particles[i].x==particles[j].x && particles[i].y==particles[j].y && particles[i].z==particles[j].z)
                    this->rep=1;
                    break;
            }
        }
    }

    else{std::cerr << "wrong data for test1";}
    this->show_rep();
}

void Test::exe(double (*Force)(particle, particle,double)){

    if (this->type==2){

        double L=10.;

        double dl=0.01;

        particle p_1(1,1,1,0,0,0);
        p_1.mass=1;
        particle p_2(9,1,1,0,0,0);
        p_2.mass=1;

        double force;
        double dist;

        std::vector<particle> particles;
        particles.push_back(p_1);
        particles.push_back(p_2);
        dist = particles[1].x - particles[0].x;

        std::ofstream force_output;
        force_output.open("./force.txt");

        while (dist>1) {
            force = (*Force)(particles[0],particles[1],L);
            dist = particles[1].x - particles[0].x;
            force_output << dist << " " <<  force << std::endl;
            particles[1].x-=dl;
        }

        force_output.close();

        system("python3 force_test_visualizer.py");
    }
    else{std::cerr << "wrong data for test2";}

}

void Test::exe(){
  if (this->type == 4){
    std::vector<double> E_vect;
    std::ofstream pos_output;
    std::ofstream energy_output;
    particle p_1(1,1,1,0,0,0);
    p_1.mass=1;
    particle p_2(3,1,1,-0.01,0,0);
    p_2.mass=1;
    std::vector<particle> Particles;
    Particles.push_back(p_1);
    Particles.push_back(p_2);
    params my_params;
    my_params.NB_PARTICLES = 2;
    my_params.NB_ITER = 500;
    pos_output.open("./positions.txt");
    energy_output.open("./energy.txt");
    pos_output <<my_params.L <<std::endl;
    pos_output <<my_params.NB_PARTICLES <<std::endl;
    pos_output <<my_params.NB_ITER <<std::endl;
    energy_output <<my_params.NB_ITER <<std::endl;
    for (int t=0; t<my_params.NB_ITER; t++){
      for (int i=0; i<my_params.NB_PARTICLES;i++){
        particle p = Particles[i];
        pos_output <<p.x <<" " <<p.y <<" " <<p.z <<std::endl;
      }
      double current_E = E(Particles, my_params.L);
      E_vect.push_back(current_E);
      energy_output  <<current_E <<std::endl;
      Stormer_Verlet(Particles, my_params.SCHEME_DT, my_params.L);
    }
    pos_output.close();
    energy_output <<mean(E_vect) <<std::endl;
    energy_output <<std_dev(E_vect) <<std::endl;
    energy_output.close();
    system("python3 visualizer.py");
  }
}
