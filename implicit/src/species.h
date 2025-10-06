#include <string>
#include <cmath>
#include <vector>
#include <functional>
#include <random>

#ifndef SPECIES_H
#define SPECIES_H

class Species
{
public:

  // numerical scheme parameters
  int m_np;
  int m_nx;
  double m_lx;
  double m_dt;
  double m_dx;

  // particle quantities
  std::string m_name;
  double m_charge;
  double m_mass;
  double m_weight;
  double* m_pos;                // np
  double* m_vel;                // np
  double* m_pos_new;            // np
  double* m_vel_new;            // np

  // moment quantities
  double* m_dens;               // row-major | ((nx) ; 2) -> dens for each cell centre for (k, k+1) timesteps
  double* m_mom;                // row-major | ((nx + 1) ; 2) -> mom for each cell face for (k+1/2) timesteps
  double* m_avgmom;             // row-major | ((nx + 1) ; 1) -> avgmom for each cell face for (k+1/2) timesteps
  double* m_stress;             // row-major | ((nx) ; 2) -> stress for each cell centre for (k, k+1) timesteps
  double* m_nstress;            // row-major | ((nx) ; 1) -> normalised stress for each cell centre for k+1/2 timestep
  double* m_gamma;              // row-major | (nx ; 1) -> consistency parameter for each cell RH face

  // LO system quantities
  double* m_elec;               // row-major | (nx ; 2)

  // ====================================================

  // tolerances
  double m_pic_tol = pow(10.0,-5);
  int m_fixed_iter_max = 4;
  int m_max_picard = 5;

public:
// empty ctor
Species();

// initialisation ctor
Species(int np, int nx, int nv, double lx, double dt, double *elec, std::string name, double charge, double mass, double v_range, double init_dens,
    std::function<double(double, double, double, double)> unscaled_pdf, std::mt19937 &generator);

// dtor
~Species();

void accumulate_moments(const double& dt_target);

// push functions
void push(const double& dt_new);

void step();

// getters
double* get_dens_ptr();
double* get_mom_ptr();
double* get_avgmom_ptr();
double* get_nstress_ptr();
double* get_gamma_ptr();

void clear_files();
void output_data();

};

#endif
