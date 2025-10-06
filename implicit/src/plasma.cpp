#include <iostream>
#include <cmath>
#include <fstream>

#include "species.h"
#include "plasma.h"

#define F77NAME(x) x##_
extern "C" {
    // LAPACK routine for solving systems of linear equations
    void F77NAME(dgesv)(const int& n, const int& nrhs, const double* A,
        const int& lda, int* ipiv, double* B,
        const int& ldb, int& info);

    // BLAS routine for performing matrix-vector multiplication
    void F77NAME(dgemv) (const char& trans, const int& m,
        const int& n, const double& alpha,
        const double* A, const int& lda,
        const double* x, const int& incx,
        const double& beta, double* y,
        const int& incy);
}

// empty ctor
Plasma::Plasma()
{
    std::cout << "Empty plasma constructor called" << std::endl;
}

// initialisation ctor
Plasma::Plasma(
    const int& _nspecies,
    const int& _nx,
    const double& _lx
)
{

    // numerical scheme parameters
    alfa_tmp = 0;
    nspecies = _nspecies;
    nx = _nx;
    lx = _lx;
    dx = lx / nx;

    // calculate x_face and x_centre;
    x_face = new double[nx + 1];
    x_centre = new double[nx];
    x_face[0] = 0.0;
    for (int i = 0; i < nx; i++)
    {
        x_face[i + 1] = (i + 1)*dx;
        x_centre[i] = 0.5 * (x_face[i] + x_face[i + 1]);
    }

    // initialise arrays
    lo_dens   = new double[nspecies * nx * 2];
    lo_avgmom = new double[nspecies * (nx+1)];
    elec      = new double[(nx+1) * 2];

    matdim = (2*nspecies + 1) * nx;
    A = new double[matdim * matdim];
    b = new double[matdim];
    soln = new double[matdim];
    lo_residual = new double[matdim];
    outer_residual = new double[matdim];
}

//dtor
Plasma::~Plasma()
{

}

// add species to plasma
void Plasma::add_species(int np, int nx, int nv, double lx, double dt, std::string name, double charge, double mass, double v_range, double init_dens,
    std::function<double(double, double, double, double)> unscaled_pdf, std::mt19937 &generator)
{
    if (alfa_tmp < nspecies)
    {
        // initialise species
        Species* species_ptr = new Species(np, nx, nv, lx, dt, this->elec, name, charge, mass, v_range, init_dens, unscaled_pdf, generator);
        
        // copy name, charge and mass to the species arrays in the Plas object
        species_name.push_back(name);
        species_charge.push_back(charge);
        species_mass.push_back(mass);
        species_ptrs.push_back(species_ptr);
        // increment species counter
        alfa_tmp++;
    }
    else
    {
        std::cout << "Adding more than " << nspecies << " species to this plasma is forbidden!" << std::endl;
    }
}

// ====================================================================
/*
initialise lo system after adding all species
  - copy HO values for current timestep k = 0 into LO values (for dens and avgmom)
  - set solution values to be same as dens and avgmom for calculating first residual
  - initialise electric field using Poisson's equation
  - set future LO values to be same as current LO values (for k = 1) for dens, elec as they are needed for A matrix
  - clean output files for density and avgmom
  - print t=0 LO density and LO elec
*/
void Plasma::init_lo()
{
    bool print_init = true;

    // use ho density and momentum for first current timestep (k = 0)
    // and set this as the solution vector for purposes of residual calculation
    for (int alfa = 0; alfa < nspecies; ++alfa)
    {
        // get ptrs to ho_dens and ho_avgmom arrays from Species object
        ho_dens = species_ptrs[alfa]->get_dens_ptr();
        ho_avgmom = species_ptrs[alfa]->get_avgmom_ptr();
        for (int i = 0; i < nx; ++i)
        {
            // lo_dens is col-major (nspecies*(nx) ; 2), ho_dens is row-major ((nx) ; 2)
            lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))] = ho_dens[2*i + 0];
            // soln is col-major ((2*nspecies + 1)*nx ; 1), ho_dens is row-major ((nx) ; 2)
            soln[i + alfa*(2*nx)] = ho_dens[2*i + 0];
            // soln is col-major ((2*nspecies + 1)*nx ; 1), ho_avgmom is row-major ((nx + 1) ; 1)
            // soln domain: [1/2 -> Nx-1/2] (RH faces) ; ho_avgmom domain: [-1/2 -> Nx-1/2] (all faces)
            // -> require i+1 for ho_avgmom indexing
            soln[i + (nx) + alfa*(2*nx)] = ho_avgmom[i+1];
        }
        for (int i = 0; i < nx + 1; ++i)
        {
            // lo_avgmom is col-major (nspecies*(nx + 1) ; 1), ho_avgmom is row-major ((nx) ; 1)
            lo_avgmom[i + alfa*(nx + 1)] = ho_avgmom[i];
        }
    }

    // estimate future lo values for first future timestep
    for (int i = 0; i < nx; ++i)
    {
        for (int alfa = 0; alfa < nspecies; ++alfa)
        {
            // lo_dens is col-major (nspecies*(nx) ; 2)
            lo_dens[cind(i + alfa*(nx), 1, nspecies*(nx))] = lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))];
        }
    }

    // initialise eletric field calculation using Poisson equation
    double* charge_dens_init = new double[nx];

    // calculate charge density
    for (int i = 0; i < nx; ++i)
    {
        charge_dens_init[i] = 0.0;
        for (int alfa = 0; alfa < nspecies; ++alfa)
        {
            double* dens_init = species_ptrs[alfa]->get_dens_ptr();
            charge_dens_init[i] += species_charge[alfa] * dens_init[2*i + 0];
        }
    }

    // integration of charge density to get E (from dE/dx = charge dens / e0)
    // E(i+1/2) = E(i-1/2) + (dt * chargeDens[i] / e0)
    elec[rind(0, 0, 2)] = 0.0;
    for (int i = 0; i < nx; ++i)
    {
        elec[rind(i + 1, 0, 2)] = charge_dens_init[i]*(dx / constants::e0) + elec[rind(i, 0, 2)];
    }
    elec[rind(0, 0, 2)] = elec[rind(0, 0, 2)] + elec[rind(nx-1, 0, 2)] - elec[rind(nx, 0, 2)];

    // calculate average electric field to make the field zero at the mid domain
    double avg_elec = 0.0;
    for (int i = 0; i < nx + 1; ++i)
    {
        avg_elec += elec[rind(i, 0, 2)];
    }
    avg_elec /= (nx + 1);

    // translate electric field to be zero at the mid domain
    for (int i = 0; i < nx + 1; ++i)
    {
        elec[rind(i, 0, 2)] -= avg_elec;
    }

    delete[] charge_dens_init;

    // initialise electric field as zeros for first current timestep
    for (int i = 0; i < nx; ++i)
    {
        // elec is row-major (nx + 1; 2) -> set RH face to zero so require i+1 index
        elec[rind(i + 1, 1, 2)] = elec[rind(i + 1, 0, 2)];
        // soln is col-major ((2*nspecies + 1)*nx ; 1)
        // elec region of soln begins at i = nspecies*(2*nx) and is nx long
        soln[i + nspecies*(2*nx)] = elec[rind(i + 1, 0, 2)];
    }
    elec[rind(0, 1, 2)] = elec[rind(0, 0, 2)];

    // initialise output files
    std::ofstream elec_file("./output/lo_elec.csv", std::ofstream::out | std::ofstream::trunc);
    elec_file.close();
    for (int alfa = 0; alfa < nspecies; ++alfa)
    {
        // for each species...
        // clear density file
        std::ofstream dens_file("./output/lo_" + species_name[alfa] + "_dens.csv", std::ofstream::out | std::ofstream::trunc);
        dens_file.close();
        // clear avgmom file
        std::ofstream avgmom_file("./output/lo_" + species_name[alfa] + "_avgmom.csv", std::ofstream::out | std::ofstream::trunc);
        avgmom_file.close();
        // clear soln file
        std::ofstream soln_file("./output/lo_" + species_name[alfa] + "_soln.csv", std::ofstream::out | std::ofstream::trunc);
        avgmom_file.close();
    }

    if (print_init == true)
    {
        // print LO dens for t=0 for each species
        for (int alfa = 0; alfa < nspecies; alfa++)
        {
            std::ofstream dens_file("./output/lo_" + species_name[alfa] + "_dens.csv", std::ofstream::out | std::ofstream::app);
            for (int i = 0; i < nx; i++)
            {
                dens_file << lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))];
                if (i < nx - 1) {
                    dens_file << ",";
                }
            }
            dens_file << "\n";
            dens_file.close();
        }

        // print electric field for t=0
        std::ofstream elec_file("./output/lo_elec.csv", std::ofstream::out | std::ofstream::app);
        for (int i = 0; i < nx + 1; i++)
        {
            elec_file << elec[rind(i, 0, 2)];
            if (i < nx) {
                elec_file << ",";
            }
        }
        elec_file << "\n";
        elec_file.close();
    }
}

// ====================================================================
/* Evolve the Plasma system using the Implicit Moment Method
  - Outer Picard loop:
    - Accumulate moments
    - Solve LO system
    - Push particles in HO systems
*/
double Plasma::evolve(const double& dt_target)
{
    // How should inconsistency error be dealt with (1 = drop time-step; 2 = recalculate electric field with Poisson with no drop in step)
    int INCONSISTENCY_FIX_METHOD = 1;

    int evolve_index = -1;

    bool evolved = false;
    int dt_reduce_factor = -1;
    double ho_lo_diff_tol = 0.01;

    // Iterative evolve loop
    while (evolved == false)
    {
        evolve_index++;
        dt_reduce_factor++;

        if (evolve_index > 5) {
            std::cout << "Stalled program as HO/LO systems won't converge" << std::endl;
            while (true) {
                
            }
        }

        // calculate timestep and set class object dt variable
        if (INCONSISTENCY_FIX_METHOD == 1)
        {
            dt = dt_target / pow(2, dt_reduce_factor);
        }
        else
        {
            dt = dt_target;
        }

        std::cout << "=================== EVOLVE ITERATION " << evolve_index << std::endl;
        std::cout << "Evolving with dt = " << dt << std::endl;

        // Set future LO values to be same as current LO values
        for (int i = 0; i < nx; i++)
        {
            for (int alfa = 0; alfa < nspecies; alfa++)
            {
                lo_dens[cind(i + alfa*(nx), 1, nspecies*(nx))] = lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))];
            }
            // requires: elec @ RH face
            elec[rind(i + 1, 1, 2)] = elec[rind(i + 1, 0, 2)];
        }
        elec[rind(0, 1, 2)] = elec[rind(0, 0, 2)];

        // initialise outer loop flag and index
        outer_flag = false;
        outer_index = -1;
        bool quit_evolve = false;

        // outer Picard loop
        while (outer_flag == false)
        {
            // increment outer index
            outer_index++;

            std::cout << "======>> Outer Picard iteration " << outer_index << std::endl;

            // =========================================
            // ACCUMULATE HO MOMENTS
            for (int alfa = 0; alfa < nspecies; alfa++)
            {
                // check target device and run corresponding function
                species_ptrs[alfa]->accumulate_moments(dt);
            }

            std::cout << "Finished momentum accumulation" << std::endl;

            // =========================================
            // SOLVE LO SYSTEM
            // initialise inner LO Picard loop flag and index
            if (outer_index > 0)
            {
                lo_flag = false;
                lo_index = -1;

                // solve LO system using inner LO Picard loop
                while (lo_flag == false)
                {
                    // increment LO index
                    lo_index++;
                    std::cout << "==>> LO inner Picard iteration " << lo_index << ": ";

                    // Update A and b matrices
                    update_A();
                    update_b();
                    
                    // solve LO system and return residual norm prior to solving to track convergence
                    lo_residual_norm = solve_lo();

                    // outer residual norm is the first residual for solving the LO system, it shows how converged the outer system is
                    if (lo_index == 0)
                    {
                        outer_residual_norm = lo_residual_norm;
                        //std::cout << "Outer residual norm: " << outer_index << ": " << outer_residual_norm << std::endl;
                    }

                    // if LO residual norm is below tolerance, terminate the inner LO Picard loop
                    if (lo_residual_norm < tol_lo_residual_norm)
                    {
                        lo_flag = true;
                    }
                    // if LO fails to converge then there is consistency error - therefore the timestep will be reduced and current evolve should be quit
                    if (lo_index > inner_max)
                    {
                        quit_evolve = true;
                        lo_flag = true;
                    }
                }
                std::cout << "LO system has been solved" << std::endl;
            }

            /*
            for (int i = 0; i < nx + 1; ++i)
            {
                std::cout << "Elec[" << i << "]: " << elec[2*i + 0] << "; " << elec[2*i + 1] << std::endl;
            }
            */

            // =========================================
            // push particles in HO system and calculate HO avgmom, unless evolve has been invalidated
            if (quit_evolve == false)
            {
                for (int alfa = 0; alfa < nspecies; alfa++)
                {
                    species_ptrs[alfa]->push(dt);
                }
            }

            std::cout << "Particles have been pushed" << std::endl;

            // =========================================
            // if the outer system has failed to converge then quit evolve as there is consistency error
            if (outer_index >= outer_max)
            {
                quit_evolve = true;
            }
            
            // check if evolve needs to be quitted
            if (quit_evolve == true)
            {
                // exit outer system and increase reduce factor
                outer_flag = true;
            }
            // check if outer system converged
            if ((outer_residual_norm < tol_outer_residual_norm) && (outer_index > 0))
            {
                outer_flag = true;
            }
        }

        // Check if LO and HO systems have converged to determine if evolve should finish
        double dens_diff = 0.0;
        double avgmom_diff = 0.0;
        for (int alfa = 0; alfa < nspecies; ++alfa)
        {
            ho_dens = species_ptrs[alfa]->get_dens_ptr();
            ho_avgmom = species_ptrs[alfa]->get_avgmom_ptr();

            for (int i = 0; i < nx; ++i)
            {
                dens_diff += fabs(lo_dens[cind(i + alfa*(nx), 1, nspecies*(nx))] - ho_dens[2*i + 1]);
                avgmom_diff += fabs(lo_avgmom[i+1 + alfa*(nx+1)] - ho_avgmom[i+1]);
            }
            avgmom_diff += fabs(lo_avgmom[0] - ho_avgmom[0]);
        }
        std::cout << "HO/LO DENS DIFF: " << dens_diff << std::endl;
        std::cout << "HO/LO AVGMOM DIFF: " << avgmom_diff << std::endl;

        // If LO and HO consistent then exit evolve loop
        if (((dens_diff + avgmom_diff)/(2*nspecies*(nx+1)) < ho_lo_diff_tol) && (quit_evolve == false))
        {
            std::cout << "Evolve completed" << std::endl;
            evolved = true;
        }
        else
        {
            if (INCONSISTENCY_FIX_METHOD == 2)
            {
                std::cout << "Recalculating electric field" << std::endl;
                // Reset problem with new eletric field recalculated with charge distribution
                // initialise eletric field calculation using Poisson equation
                double* charge_dens_init = new double[nx];

                // calculate charge density
                for (int i = 0; i < nx; ++i)
                {
                    charge_dens_init[i] = 0.0;
                    for (int alfa = 0; alfa < nspecies; ++alfa)
                    {
                    double* dens_init = species_ptrs[alfa]->get_dens_ptr();
                    charge_dens_init[i] += species_charge[alfa] * dens_init[2*i + 0];
                    }
                }

                // integration of charge density to get E (from dE/dx = charge dens / e0)
                // E(i+1/2) = E(i-1/2) + (dt * chargeDens[i] / e0)
                elec[rind(0, 0, 2)] = 0.0;
                for (int i = 0; i < nx; ++i)
                {
                    elec[rind(i + 1, 0, 2)] = charge_dens_init[i]*(dx / constants::e0) + elec[rind(i, 0, 2)];
                }
                elec[rind(0, 0, 2)] = elec[rind(0, 0, 2)] + elec[rind(nx-1, 0, 2)] - elec[rind(nx, 0, 2)];

                // calculate average electric field to make the field zero at the mid domain
                double avg_elec = 0.0;
                for (int i = 0; i < nx + 1; ++i)
                {
                    avg_elec += elec[rind(i, 0, 2)];
                }
                avg_elec /= (nx + 1);

                // translate electric field to be zero at the mid domain
                for (int i = 0; i < nx + 1; ++i)
                {
                    elec[rind(i, 0, 2)] -= avg_elec;
                }

                delete[] charge_dens_init;
            }
        }
    }

    return dt;
}

// ====================================================================
void Plasma::update_A()
{
    //reset A matrix to zeros
    for (int i = 0; i < matdim*matdim; i++)
    {
        A[i] = 0.0;
    }

    // Loop through species for submatrices for species-dependent quantities
    for (int alfa = 0; alfa < nspecies; alfa++)
    {
        // get required HO quantities (all arranged row-major)
        ho_nstress = species_ptrs[alfa]->get_nstress_ptr();
        ho_gamma = species_ptrs[alfa]->get_gamma_ptr();

        // loop through cell indexes
        for (int i = 0; i < nx; i++)
        {
        // coefficients that do depend on species
        // density for eq. 2
        if (i < nx - 1)
        {
            // requires: dens and nstress for current + future timesteps (k, k+1) for this and RH cells
            A[cind(i + (nx) + alfa*(2*nx), i + alfa*(2*nx), matdim)] = 0.5*(-(ho_nstress[i]/dx) - ho_gamma[i]);
            A[cind(i + (nx) + alfa*(2*nx), i+1 + alfa*(2*nx), matdim)] = 0.5*((ho_nstress[i+1]/dx) - ho_gamma[i]);
        }
        else
        {
            // for cell i=nx-1, periodic boundary conditions: nstress from RH cell must be from cell i=0 instead
            A[cind(i + (nx) + alfa*(2*nx), i + alfa*(2*nx), matdim)] = 0.5*(-(ho_nstress[i]/dx) - ho_gamma[i]);
            A[cind(i + (nx) + alfa*(2*nx), 0 + alfa*(2*nx), matdim)] = 0.5*((ho_nstress[0]/dx) - ho_gamma[i]);
        }

        // electric field for eq. 2
        if (i < nx - 1)
        {
            // requires: dens for current + future timesteps (k, k+1) for this and RH cells
            A[cind(i + (nx) + alfa*(2*nx), i + nspecies*(2*nx), matdim)] = lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))] + lo_dens[cind(i+1 + alfa*(nx), 0, nspecies*(nx))]
            + lo_dens[cind(i + alfa*(nx), 1, nspecies*(nx))] + lo_dens[cind(i+1 + alfa*(nx), 1, nspecies*(nx))];
            A[cind(i + (nx) + alfa*(2*nx), i + nspecies*(2*nx), matdim)] *= -0.125 * (species_charge[alfa] / species_mass[alfa]);
        }
        else
        {
            // for cell i=nx-1, periodic boundary conditions: dens from RH cell must be from cell i=0 instead
            A[cind(i + (nx) + alfa*(2*nx), i + nspecies*(2*nx), matdim)] = lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))] + lo_dens[cind((0) + alfa*(nx), 0, nspecies*(nx))]
            + lo_dens[cind(i + alfa*(nx), 1, nspecies*(nx))] + lo_dens[cind((0) + alfa*(nx), 1, nspecies*(nx))];
            A[cind(i + (nx) + alfa*(2*nx), i + nspecies*(2*nx), matdim)] *= -0.125 * (species_charge[alfa] / species_mass[alfa]);
        }
        // avg momentum for eq. 3
        A[cind(i + nspecies*(2*nx), i + (nx) + alfa*(2*nx), matdim)] = species_charge[alfa];

        // coefficients that do not depend on species
        // density for eq. 1
        A[cind(i + alfa*(2*nx), i + alfa*(2*nx), matdim)] = 1.0 / dt;

        // avg momentum for eq. 1
        A[cind(i + alfa*(2*nx), i + (nx) + alfa*(2*nx), matdim)] = 1.0 / dx;
        if (i > 0)
        {
            A[cind(i + alfa*(2*nx), i-1 + (nx) + alfa*(2*nx), matdim)] = -1.0 / dx;
        } else {
            A[cind(i + alfa*(2*nx), nx-1 + (nx) + alfa*(2*nx), matdim)] = -1.0 / dx;
        }

        // avg momentum for eq. 2
        A[cind(i + (nx) + alfa*(2*nx), i + (nx) + alfa*(2*nx), matdim)] = 2.0 / dt;
        }
    }

    for (int i = 0; i < nx; i++)
    {
        // electric field for eq. 3
        A[cind(i + nspecies*(2*nx), i + nspecies*(2*nx), matdim)] = constants::e0 / dt;
    }

    // periodic boundary condition for electric field gradient

    // print A matrix
/*
  std::ofstream A_opfile("./output/A.csv", std::ofstream::out | std::ofstream::app);
  for (int i = 0; i < matdim; i++) {
    for (int j = 0; j < matdim; j++) {
      A_opfile.width(15);
      A_opfile << A[cind(i, j, matdim)] << "\t";
    }
    A_opfile << std::endl;
  }
  A_opfile << std::endl;
  A_opfile << "########################################################################" << std::endl;
  A_opfile << std::endl;
  A_opfile.close();
 */
}

// ====================================================================
// update b vector
void Plasma::update_b() {
    for (int alfa = 0; alfa < nspecies; alfa++)
    {
        // get pointers to required HO moments
        ho_dens = species_ptrs[alfa]->get_dens_ptr();
        ho_mom = species_ptrs[alfa]->get_mom_ptr();
        ho_nstress = species_ptrs[alfa]->get_nstress_ptr();

        for (int i = 0; i < nx; i++)
        {
        // eq. 1
        b[i + alfa*(2*nx)] = ho_dens[2*i + 0] / dt;

        // eq. 2
        if (i < nx - 1)
        {
            // requires: mom @ RH face for each cell
            b[i + (nx) + alfa*(2*nx)] = 2.0 * ho_mom[2*(i+1) + 0] / dt;
            // requires: dens for current timestep (k) and nstress from this and RH cells
            b[i + (nx) + alfa*(2*nx)] += (0.5 / dx) * ((lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))] * ho_nstress[i]) - (lo_dens[cind(i+1 + alfa*(nx), 0, nspecies*(nx))] * ho_nstress[i+1]));
            // requires: dens for current + future timestep (k, k + 1) from this and RH cells
            // requires: elec for current timestep (k) from RH face, use i+1 index
            b[i + (nx) + alfa*(2*nx)] += 0.125 * (species_charge[alfa] / species_mass[alfa]) * (lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))] + lo_dens[cind(i+1 + alfa*(nx), 0, nspecies*(nx))]
                + lo_dens[cind(i + alfa*(nx), 1, nspecies*(nx))] + lo_dens[cind(i+1 + alfa*(nx), 1, nspecies*(nx))]) * elec[rind(i + 1, 0, 2)];
        }
        else
        {
            // mom @ RH face for each cell
            b[i + (nx) + alfa*(2*nx)] = 2.0 * ho_mom[2*(i+1) + 0] / dt;
            // requires: dens and nstress from RH cell @ i=Nx-1 must come from cell i=0
            b[i + (nx) + alfa*(2*nx)] += (0.5 / dx) * ((lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))] * ho_nstress[i]) - (lo_dens[cind((0) + alfa*(nx), 0, nspecies*(nx))] * ho_nstress[(0)]));
            // requires: dens from RH cell @ i=Nx-1 must come from cell i=0
            // requires: elec for current timestep (k) from RH face, use i+1 index
            b[i + (nx) + alfa*(2*nx)] += 0.125 * (species_charge[alfa] / species_mass[alfa]) * (lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))] + lo_dens[cind((0) + alfa*(nx), 0, nspecies*(nx))]
                + lo_dens[cind(i + alfa*(nx), 1, nspecies*(nx))] + lo_dens[cind((0) + alfa*(nx), 1, nspecies*(nx))]) * elec[rind(i + 1, 0, 2)];
        }
        }
    }
    for (int i = 0; i < nx; i++)
    {
        // eq. 3
        // requires: elec @ RH face for each cell - use i+1 index
        b[i + nspecies*(2*nx)] = constants::e0 * elec[rind(i + 1, 0, 2)] / dt;

        // NEW : averaged component (using previous Picard iteration)
        for (int alfa = 0; alfa < nspecies; ++alfa)
        {
            for (int j = 0; j < nx + 1; ++j)
            {
                b[i + nspecies*(2*nx)] += species_charge[alfa] * (1.0 / (nx + 1)) * (lo_avgmom[j + alfa*(nx+1)]);
            }
        }
    }


    // print b file
    /*
    std::ofstream b_opfile("./output/b.csv", std::ofstream::out | std::ofstream::app);
    for (int i = 0; i < matdim; i++) {
        b_opfile.width(15);
        b_opfile << b[i];
        b_opfile << std::endl;
    }
    b_opfile << std::endl;
    b_opfile << "########################################################################" << std::endl;
    b_opfile << std::endl;
    b_opfile.close();
*/
}

// ====================================================================
// solve LO system and return convergence
double Plasma::solve_lo()
{
    // Calculate LO residual using prevous solution vector (either previous LO iteration or previous outer iteration)
    // transfer data from b vector to soln vector
    for (int i = 0; i < matdim; ++i)
    {
        lo_residual[i] = b[i];
    }

    // perform matrix-vector multiplication and addition to b vector
    //cblas_dgemv(CblasColMajor, CblasNoTrans, matdim, matdim, 1.0, A, matdim, soln, 1, -1.0, lo_residual, 1);
    F77NAME(dgemv)('N', matdim, matdim, 1.0, A, matdim, soln, 1, -1.0, lo_residual, 1);

    // calculate residual norm
    double residual_norm = 0.0;
    for (int i = 0; i < matdim; ++i)
    {
        residual_norm += fabs(lo_residual[i]);
    }
    std::cout << "Residual norm: " << residual_norm << std::endl;

    // initialise solve variables
    int info = 0;
    int nrhs = 1;
    ipiv = new int[matdim];

    // copy b values to soln vector
    for (int i = 0; i < matdim; i++)
    {
        soln[i] = b[i];
    }

    // solve for new solution using LAPACK library dgesv function
    F77NAME(dgesv)(matdim, nrhs, A, matdim, ipiv, soln, matdim, info);

    // update future values
    for (int i = 0; i < nx; i++)
    {
        for (int alfa = 0; alfa < nspecies; alfa++)
        {
        // requires: dens for future timestep (k+1) for each cell
        lo_dens[cind(i + alfa*(nx), 1, nspecies*(nx))] = soln[i + alfa*(2*nx)];
        // requires: avgmom for mid timestep (k+1/2) for RH face of each cell -> use i+1 index
        lo_avgmom[i + 1 + alfa*(nx+1)] = soln[i + (nx) + alfa*(2*nx)];
        }
        // requires: elec for future timestep (k+1) for RH face of each cell -> use i+1 index
        elec[rind(i + 1, 1, 2)] = soln[i + nspecies*(2*nx)];
    }

    // periodic boundary condition for LO avgmom and electric field -> LH face of cell i=0 = RH face of cell i=nx-1
    for (int alfa = 0; alfa < nspecies; alfa++)
    {
        // requires: avgmom for mid timestep (k+1/2) for LH face of cell i=0
        lo_avgmom[0 + alfa*(nx+1)] = lo_avgmom[nx + alfa*(nx+1)];
    }
    // requires: elec for future timestep (k+1) for LH face of cell i=0
    // gradient boundary condition at edges - gradient of electric field (potential) should be equal
    elec[rind(0, 1, 2)] = -elec[rind(nx, 1, 2)] + elec[rind(nx-1, 1, 2)] + elec[rind(1, 1, 2)];

    // warn if there is a singular matrix
    if (info != 0)
    {
        std::cout << "Singular A matrix in LO system!" << std::endl;
    }

    delete[] ipiv;

    //print_vals(0, 1, true);

    return residual_norm;
}

// ====================================================================
// forward step LO system
void Plasma::step() {
    // equate density for old future timestep to HO density for old future timestep and move to new current timestep
    // set solution vector to be same as current values for residual
    for (int alfa = 0; alfa < nspecies; alfa++)
    {
        ho_dens = species_ptrs[alfa]->get_dens_ptr();
        ho_avgmom = species_ptrs[alfa]->get_avgmom_ptr();
        for (int i = 0; i < nx; i++) {
            // requires: dens for current timestep (k) (soon to be future timestep) for each cell centre
            lo_dens[cind(i + alfa*(nx), 0, nspecies*(nx))] = ho_dens[2*i + 1];
            // requires: avgmom for mid timestep (k+1/2) -> Is this even necessary?
            lo_avgmom[i+1 + alfa*(nx+1)] = ho_avgmom[i+1];

            // use old future values for future solution vector for calculating initial outer residual
            soln[i + alfa*(2*nx)] = ho_dens[2*i + 1];
            soln[i + (nx) + alfa*(2*nx)] = ho_avgmom[i+1];
        }
        lo_avgmom[0 + alfa*(nx+1)] = ho_avgmom[0];
    }

    // transfer electric field from old future timestep to new current timestep
    for (int i = 0; i < nx; i++)
    {
        elec[rind(i + 1, 0, 2)] = elec[rind(i + 1, 1, 2)];

        // use old future electric field for future solution vector for calculating initial outer residual
        soln[i + nspecies*(2*nx)] = elec[rind(i + 1, 1, 2)];
    }
    // periodic boundary condition: elec @ LH face of cell i=0
    elec[rind(0, 0, 2)] = elec[rind(0, 1, 2)];

    // forward step species
    for (int alfa = 0; alfa < nspecies; alfa++)
    {
        species_ptrs[alfa]->step();
    }
}

// ====================================================================
// indexing functions
int Plasma::rind(const int& row, const int& col, const int& ncols)
{
    return row*ncols + col;
}

int Plasma::cind(const int& row, const int& col, const int& nrows)
{
    return row + col*nrows;
}

// output function
void Plasma::output_data()
{
    // print LO electric field data
    std::ofstream elec_file("./output/lo_elec.csv", std::ofstream::out | std::ofstream::app);
    for (int i = 0; i < nx + 1; i++)
    {
        elec_file << elec[rind(i, 1, 2)];
        if (i < nx) {
            elec_file << ",";
        }
    }
    elec_file << "\n";
    elec_file.close();

    for (int alfa = 0; alfa < nspecies; alfa++)
    {
        // print LO density data
        std::ofstream dens_file("./output/lo_" + species_name[alfa] + "_dens.csv", std::ofstream::out | std::ofstream::app);
        for (int i = 0; i < nx; i++)
        {
            dens_file << lo_dens[cind(i + alfa*(nx), 1, nspecies*(nx))];
            if (i < nx - 1) {
                dens_file << ",";
            }
        }
        dens_file << "\n";
        dens_file.close();

        // print LO avgmom data
        std::ofstream avgmom_file("./output/lo_" + species_name[alfa] + "_avgmom.csv", std::ofstream::out | std::ofstream::app);
        for (int i = 0; i < nx + 1; i++)
        {
            avgmom_file << lo_avgmom[i + alfa*(nx+1)];
            if (i < nx) {
                avgmom_file << ",";
            }
        }
        avgmom_file << "\n";
        avgmom_file.close();
    }

    // print LO system solution vector
    for (int alfa = 0; alfa < nspecies; alfa++)
    {
        std::ofstream soln_file("./output/lo_" + species_name[alfa] + "_soln.csv", std::ofstream::out | std::ofstream::app);
        for (int i = 0; i < nx; i++)
        {
            soln_file << soln[i + alfa*(2*nx)];
            if (i < nx - 1) {
                soln_file << ",";
            }
        }
        for (int i = 0; i < nx; i++)
        {
            soln_file << soln[i + (nx) + alfa*(2*nx)];
            if (i < nx - 1) {
                soln_file << ",";
            }
        }
        soln_file << "\n";
        soln_file.close();
    }

    // print values from species
    for (int alfa = 0; alfa < nspecies; alfa++)
    {
        species_ptrs[alfa]->output_data();
    }

    std::cout << "Species printed values" << std::endl;
}