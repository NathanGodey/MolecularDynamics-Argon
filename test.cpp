#include "params.h"
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
        particle p_2(1+dl*500,1,1,-dl,0,0);
        p_2.mass=1;

        double force;
        double dist;

        std::vector<particle> particles;
        particles.push_back(p_1);
        particles.push_back(p_2);
        dist = particles[1].x - particles[0].x;

        std::ofstream force_output;
        force_output.open("/home/adrien/COURS_ENPC/deuxieme_annee/Projet_Phy_code/force.txt");

        while (dist>0.01) {
            force = (*Force)(particles[0],particles[1],L);
            dist = particles[1].x - particles[0].x;
            force_output << dist << " " <<  force << std::endl;
            particles[1].evolve(L,0.01);
        }

        force_output << -1 << " " << 0 << std::endl;
        force_output.close();

        system("python3 /home/adrien/COURS_ENPC/deuxieme_annee/Projet_Phy_code/visualizer_test2.py");
    }
    else{std::cerr << "wrong data for test2";}

}
