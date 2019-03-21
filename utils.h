#pragma once
#include <vector>


class particle {
public:
  double x, y, z, px, py, pz;
  double mass;
  particle();
  particle(double x, double y, double z, double px, double py, double pz);
  void evolve(double , double dt);
};
