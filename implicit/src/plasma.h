#include <string>
#include <cmath>
#include <vector>

#include "species.h"

namespace constants
{
    const double e0 = 1.0;//8.85418 * pow(10.0, -12);
	const double Kb = 1.0;//1.3806 * pow(10.0, -23);
}

#ifndef PLASMA_H
#define PLASMA_H

class Plasma
{
public:

  //numerical scheme parameters
  int nspecies;
  int nx;
  double lx;
  double dx;
  double dt;
  double* x_centre;
  double* x_face;

  // species quantities
  std::vector<std::string> species_name;
  std::vector<double> species_charge;
  std::vector<double> species_mass;
  std::vector<Species*> species_ptrs;
  int alfa_tmp;

  // LO quantities
  double* lo_dens;                // col-major | (nx*(nspecies + 1) ; 2) -> LO density for each cell centre for current + future timestep (k, k+1)
  double* lo_avgmom;              // col-major | ((nx + 1)*(nspecies) ; 1) -> LO averge momentum for each cell face for mid timestep (k+1/2)
  double* elec;                   // row-major | ((nx + 1) ; 2) -> electric field for each cell face for current + future timestep (k, k+1)

  // HO quantities
  double* ho_dens;                // row-major | ((nx) ; 2) -> HO density for each cell centre for current + future timestep (k, k+1)
  double* ho_mom;                 // row-major | ((nx + 1) ; 1) -> HO momentum for each cell face for current timestep (k)
  double* ho_avgmom;              // row-major | ((nx + 1) ; 1) -> HO average momentum for each cell face for mid timestep (k+1/2)
  double* ho_nstress;             // row-major | ((nx) ; 1) -> HO normalised stress for each cell centre for mid timestep (k+1/2)
  double* ho_gamma;               // row-major | ((nx) ; 1) -> consistency parameter for RH face of each cell for mid timestep (k+1/2)

  // calculation variables
  double* A;                      // col-major | (nx*(2*nspecies + 1) ; nx*(2*nspecies + 1)) -> A matrix for solving LO
  double* b;                      // col-major | (nx*(2*nspecies + 1) ; 1) -> b vector for solving LO
  double* soln;                   // col-major | (nx*(2*nspecies + 1) ;'1) -> solution vector to contain LO solution
  double* outer_residual;         // col-major | (nx*(2*nspecies + 1) ; 1) -> residual vector
  double outer_residual_norm;
  double* lo_residual;
  double lo_residual_norm;
  double tol_lo_residual_norm     = pow(10.0, -6);
  double tol_outer_residual_norm  = pow(10.0, -6);
  int outer_max = 30;
  int inner_max = 30;
  int matdim;
  int* ipiv;

  bool lo_flag;
  bool outer_flag;
  int lo_index;
  int outer_index;

public:
  Plasma();
  Plasma(
    const int& _nspecies,
    const int& _nx,
    const double& _lx
  );

  ~Plasma();

  void add_species(int np, int nx, int nv, double lx, double dt, std::string name, double charge, double mass, double v_range, double init_dens,
    std::function<double(double, double, double, double)> unscaled_pdf, std::mt19937 &generator);
  void init_lo();

  double evolve(const double& dt_target);
  void update_A();
  void update_b();
  double solve_lo();
  void step();
  int rind(const int& row, const int& col, const int& ncols);
  int cind(const int& row, const int& col, const int& nrows);
  void output_data();
};

#endif
