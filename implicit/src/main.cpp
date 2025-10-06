#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <random>

#include "species.h"
#include "plasma.h"
//#include "config.cpp"

double unscaled_pdf_landau(double length, double v_range, double x, double v) {
    double fx = 1 + (0.2 * sin(2.0 * M_PI * x / length));
    double fv = std::exp(-v*v / (2 * (v_range*v_range / 1.96*1.96)));
    return fx * fv;
}

double uniform_pdf(double length, double v_range, double x, double v) {
  return 1.0;
}

int main()
{

  //==================================================================
  // config objects for different benchmark tests

/*
  Config landau1;
  landau1.nspecies = 1;
  landau1.np = (int) 10000;
  landau1.nx = 30;
  landau1.nv = { 50 , 50 };
  landau1.lx = 4.0 * M_PI;
  landau1.dt = 1.0 * 0.01 * 1.467;
  landau1.time_total = 200;
  landau1.species_name = { "electron", "ion" };
  landau1.species_dens_avg = { 1.0, 1.0 };
  landau1.species_dens_perturb = { 0.1, 0.0 };
  landau1.species_dens_perturb_profile = { "cos", "" };
  landau1.species_vel_profile = { "boltzmann", "two_beams" };
  landau1.species_vel_avg = { 0.0, 0.0 };
  landau1.species_vel_range = { 8.0 , 0.0 };
  landau1.species_vel_perturb = { 0.0, 0.0 };
  landau1.species_vel_perturb_profile = { "", "" };
  landau1.species_charge = { -1.0 , 1.0 };
  landau1.species_mass = { 1.0 , 2000.0 };
  landau1.species_T = { 1.0 , 1.0 };
*/
  //==================================================================

/*
  int skip = 1;
  int test_config_index = 1;

  // use direct method/fixed method accelerator?
  int accelerate = 0;

  // create plasma object and associated species objects
  landau1.create_plasma();

  std::vector<double>runtime_data = landau1.run_plasma(accelerate, skip);

  std::cout << "PROGRAM COMPLETE" << std::endl;
  std::cout << "Executed in " << runtime_data[0] << " seconds" << std::endl;
  std::cout << "Simulation time: " << runtime_data[1] << " seconds" << std::endl;

*/

// seed 42 causes issues
  int seed = 43;

  std::random_device rd;
  std::mt19937 generator(seed);

  int nx = 16;
  double lx = 4.0 * M_PI;
  int np = 2 << 16;
  double dt = 0.01 * 1.467;
  double electron_charge = -1.0;
  double electron_mass = 1.0;
  double ion_charge = 1.0;
  double ion_mass = 2000.0;

  Plasma *plasma;

  plasma = new Plasma(2, nx, lx);

  plasma->add_species(np, nx, 8, lx, dt, "electron", electron_charge, electron_mass, 8.0, 1.0, unscaled_pdf_landau, generator);
  plasma->add_species(np, nx, 2, lx, dt, "ion", ion_charge, ion_mass, 0.001, 1.0, uniform_pdf, generator);

  plasma->init_lo();

  plasma->species_ptrs[0]->push(dt);

  plasma->species_ptrs[0]->accumulate_moments(dt);

  for (int i = 0; i < nx + 1; ++i) {
    std::cout << plasma->species_ptrs[0]->m_avgmom[i] << ", " << plasma->species_ptrs[0]->m_mom[2*i + 1] << std::endl;
  }


  for (int k = 0; k < 300; ++k) {
    std::cout << "========= SOLVING TIMESTEP " << k+1 << "===========" << std::endl;
    plasma->evolve(dt);

    if (k % 5 == 0) {
      plasma->output_data();
    }

    plasma->step();
  }

  return 0;
}
