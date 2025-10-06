#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <functional>
#include <random>

#include "species.h"

// sign function (imported)
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// empty ctor
Species::Species()
{
    std::cout << "Empty species constructor called" << std::endl;
}

// helper function for particle initialisation
std::vector<double> get_pdf_grid(int Ncells, int Nvcells, double length, double v_range, std::function<double(double,double,double,double)> unscaled_pdf) {
    std::vector<double> pdf(Ncells * Nvcells);
    double dx = length / Ncells;
    double dv = 2 * v_range / Nvcells;
    double total = 0;
    
    for (int i = 0; i < Ncells; ++i) {
        for (int j = 0; j < Nvcells; ++j) {

            double xcell_center = (0.5 + i) * dx;
            double vcell_center = -v_range + (0.5 + j) * dv;

            pdf[i*Nvcells + j] = unscaled_pdf(length, v_range, xcell_center, vcell_center);
            
            total += pdf[i*Nvcells + j];
        }
    }
    for (int i = 0; i < Ncells * Nvcells; ++i) {
        pdf[i] /= total;
    }
    return pdf;
}

// initialisation ctor
Species::Species(int np, int nx, int nv, double lx, double dt, double *elec, std::string name, double charge, double mass, double v_range, double init_dens,
    std::function<double(double, double, double, double)> unscaled_pdf, std::mt19937 &generator)
{
    // copy parameters to class member variables
    m_np = np;
    m_nx = nx;
    m_lx = lx;
    m_dt = dt;
    m_elec  = elec;
    m_name = name;
    m_charge = charge;
    m_mass = mass;
    m_dx = lx / nx;

    // allocate memory for arrays
    m_pos     = new double[m_np];
    m_pos_new = new double[m_np];
    m_vel     = new double[m_np];
    m_vel_new = new double[m_np];
    m_dens    = new double[m_nx * 2];
    m_mom     = new double[(m_nx + 1) * 2]; // does this need to be for current timestep and future timestep?
    m_avgmom  = new double[m_nx + 1];
    m_stress  = new double[m_nx * 2];
    m_nstress = new double[m_nx];
    m_gamma   = new double[m_nx];

    // reset output files
    clear_files();

    // now perform particle initialisation

    double dx = m_lx / m_nx;
    double dv = 2 * v_range / nv;

    // create a uniform distribution object which be sampled from to distribute particles
    // the x-v space is discretised into cells and the pdf is used to calculate a target number
    // of particles for each x-v cell. the distribution of particles in each x-v cell is uniform
    std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

    // get the pdf grid -- that is a grid of weightings for each x-v cell that sum to 1 and indicate
    // how many particles will be created in each cell. the pdf grid has dimensions nx * nv
    // and for x=i, v=j (by indices) it is accessed via pdf_grid[i*nv + j]
    std::vector<double> pdf_grid = get_pdf_grid(m_nx, nv, m_lx, v_range, unscaled_pdf);

    int particle_idx = 0;

    for (int i = 0; i < m_nx; ++i) {
        for (int j = 0; j < nv; ++j) {
            // calculate the target number of particles for this x-v cell
            // this rounds down, therefore there is no risk of "overallocating" the particles
            int target_np = int(np * pdf_grid[i*nv + j]);
            for (int k = 0; k < target_np; ++k) {
                // particles are uniformly distributed in this x-v cell
                double x = (i + uniform_distribution(generator)) * dx;
                double v = -v_range + (j + uniform_distribution(generator)) * dv;

                m_pos[particle_idx] = x;
                m_vel[particle_idx] = v;

                particle_idx++;
            }
        }
    }

    // randomly place the remaining particles that need to be initialised
    for ( ; particle_idx < m_np; ++particle_idx) {
        // generate a random x,v value
        double x = uniform_distribution(generator) * m_lx;
        double v = -v_range + uniform_distribution(generator) * 2 * v_range;
        m_pos[particle_idx] = x;
        m_vel[particle_idx] = v;
    }

    // now calculate the current density, momentum and stress

    // set current dens, momentum and stress to zero ahead of accumulation
    for (int i = 0; i < m_nx; ++i) {
        m_dens[2*i + 0]   = 0.0;
        m_stress[2*i + 0] = 0.0;
    }
    for (int i = 0; i < m_nx + 1; ++i) {
        m_mom[2*i + 0] = 0.0;
    }

    for (int p = 0; p < m_np; ++p) {
        // get the root index for the cell that this particle is contained in
        // this index corresponds to the lhs face of the cell and the index of the cell centre
        int cell_idx = (int) floor(m_pos[p] / m_dx);

        // calculate first-order weights for current position
        double weight_lhs_face = 1.0 - (m_pos[p] - cell_idx*m_dx) / m_dx;
        double weight_rhs_face = 1.0 - weight_lhs_face;

        // calculate the second-order weights for current position
        // first calculate |x-x_p| / dx for the lhs, middle and rhs cell centre positions
        double d_lhs_centre = fabs(m_pos[p] - ((double) cell_idx - 0.5)*m_dx) / m_dx;
        double d_mid_centre = fabs(m_pos[p] - ((double) cell_idx + 0.5)*m_dx) / m_dx;
        double d_rhs_centre = fabs(m_pos[p] - ((double) cell_idx + 1.5)*m_dx) / m_dx;

        // for middle cell centre, always <= 0.5
        double weight_mid_centre = 0.75 - d_mid_centre*d_mid_centre;
        // for LH and RH neighbouring cell centres, always > 0.5
        double weight_lhs_centre = 0.5 * (1.5 - d_lhs_centre) * (1.5 - d_lhs_centre);
        double weight_rhs_centre = 0.5 * (1.5 - d_rhs_centre) * (1.5 - d_rhs_centre);

        // get the indexes for the lhs and rhs cells for the density and stress accumulation
        int lhs_cell_idx = cell_idx - 1;
        if (lhs_cell_idx < 0) {
            lhs_cell_idx = m_nx - 1;
        }
        int rhs_cell_idx = (cell_idx + 1) % m_nx;

        // accumulate moment contributions (without weight / dx contribution)
        // initialise density @ t = 0
        m_dens[2*cell_idx     + 0] += weight_mid_centre;
        m_dens[2*lhs_cell_idx + 0] += weight_lhs_centre;
        m_dens[2*rhs_cell_idx + 0] += weight_rhs_centre;

        // initialise momentum @ t = 0
        m_mom[2*cell_idx       + 0] += weight_lhs_face * m_vel[p];
        m_mom[2*(cell_idx + 1) + 0] += weight_rhs_face * m_vel[p]; // change this to be rhs_cell_idx, when reformulating momentum to be nx long instead of (nx + 1)

        // initialise stress @ t = 0
        m_stress[2*lhs_cell_idx + 0] += weight_lhs_centre * m_vel[p] * m_vel[p];
        m_stress[2*cell_idx     + 0] += weight_mid_centre * m_vel[p] * m_vel[p];
        m_stress[2*rhs_cell_idx + 0] += weight_rhs_centre * m_vel[p] * m_vel[p];
    }

    // finish calculating the density in the case that the weight = 1.0
    for (int i = 0; i < m_nx; ++i) {
        m_dens[2*i + 0] *= 1.0 / m_dx;
    }

    // calculate avg density when weight = 1.0 for each particle
    // then this will be used to calculate the weight such that the average density is equal to init_dens
    double avg_dens = 0.0;
    for (int i = 0; i < m_nx; ++i) {
        avg_dens += m_dens[2*i + 0];
    }
    avg_dens /= m_nx;

    // scale particle weight so avg density matches target density
    m_weight = init_dens / avg_dens;
    
    for (int i = 0; i < m_nx; ++i)
    {
        // now scale density by the calculated weight
        m_dens[2*i + 0] *= m_weight;

        // finish off the calculation of stress
        m_stress[2*i + 0] *= m_weight / m_dx;
    }
    for (int i = 0; i < m_nx + 1; ++i)
    {
        // finish off the calculation of momentum
        m_mom[2*i + 0] *= m_weight / m_dx;
    }

    // periodic boundary condition for momentum
    m_mom[2*0    + 0] += m_mom[2*m_nx + 0];
    m_mom[2*m_nx + 0] =  m_mom[2*0    + 0]; // this will not be needed when changing momentum to only be nx longer vector

    // set avgmom to zeros for mid timestep 1/2
    std::fill(m_avgmom, m_avgmom + m_nx + 1, 0.0);

    // copy the current position and velocity into future position and velocity
    std::copy(m_pos, m_pos + m_np, m_pos_new);
    std::copy(m_vel, m_vel + m_np, m_vel_new);
}

//==============================================================================
//==============================================================================

/* update moments:
-> dens (time k+1),
-> mom (time k),
-> stress (time k+1),
-> nstress (time k+1/2),
-> gamma (time k+1/2)
*/
void Species::accumulate_moments(const double& dt_target)
{
    m_dt = dt_target;

    // set relevant vectors to zero ahead of accumulation
    for (int i = 0; i < m_nx; ++i)
    {
        // set future density to zero
        m_dens[2*i + 1] = 0.0;

        // set future stress to zero
        m_stress[2*i + 1] = 0.0;
    }
    for (int i = 0; i < m_nx + 1; ++i)
    {
        // set current momentum to zero            [THIS IS AN ERROR? SHOULD BE 2*i + 0 for current momentum, WE ARE CALCULATING FUTURE MOMENTUM!]
        m_mom[2*i + 1] = 0.0;
    }

    //******************************************************************
    // accumulation step -- loop through particles
    for (int p = 0; p < m_np; ++p)
    {
        // calculate LH and RH cell wall indices for current position
        int cell_lhface = (int) floor(m_pos_new[p] / m_dx);
        int cell_rhface = cell_lhface + 1;

        // calculate LH, middle and RH face indices for new position
        int cell_midcentre = (int) floor(m_pos_new[p] / m_dx);
        int cell_lhcentre = cell_midcentre - 1;
        int cell_rhcentre = cell_midcentre + 1;

        // calculate first-order shapes for current position
        double shape_lhface = 1.0 - fabs(m_pos_new[p] / m_dx - (double) cell_lhface);
        double shape_rhface = 1.0 - fabs(m_pos_new[p] / m_dx - (double) cell_rhface);

        // calculate the second-order shapes for future new position
        double d_lhcentre  = fabs(m_pos_new[p] / m_dx - ((double) cell_lhcentre + 0.5));
        double d_midcentre = fabs(m_pos_new[p] / m_dx - ((double) cell_midcentre + 0.5));
        double d_rhcentre  = fabs(m_pos_new[p] / m_dx - ((double) cell_rhcentre + 0.5));

        // calculate second order shapes
        // for middle cell centre, always <= 0.5
        double shape_midcentre = 0.75 - (d_midcentre * d_midcentre);
        // for LH and RH neighbouring cell centres, always > 0.5
        double shape_lhcentre = 0.5 * (1.5 - d_lhcentre) * (1.5 - d_lhcentre);
        double shape_rhcentre = 0.5 * (1.5 - d_rhcentre) * (1.5 - d_rhcentre);

        // apply periodic boundary conditions to indices for accumulation
        // face indices
        if (cell_lhcentre < 0)
        {
            cell_lhcentre = m_nx - 1;
        }
        cell_rhcentre = cell_rhcentre % m_nx;

        // accumulate moment contributions (without weight / dx contribution)
        // density for t = k+1
        m_dens[2*cell_lhcentre  + 1] += shape_lhcentre;
        m_dens[2*cell_midcentre + 1] += shape_midcentre;
        m_dens[2*cell_rhcentre  + 1] += shape_rhcentre;

        // momentum for t = k+1
        m_mom[2*cell_lhface + 1] += m_vel_new[p] * shape_lhface;
        m_mom[2*cell_rhface + 1] += m_vel_new[p] * shape_rhface;

        // stress for t = k+1
        m_stress[2*cell_lhcentre  + 1] += m_vel_new[p] * m_vel_new[p] * shape_lhcentre;
        m_stress[2*cell_midcentre + 1] += m_vel_new[p] * m_vel_new[p] * shape_midcentre;
        m_stress[2*cell_rhcentre  + 1] += m_vel_new[p] * m_vel_new[p] * shape_rhcentre;
    }

    //******************************************************************

    for (int i = 0; i < m_nx; ++i)
    {
        m_dens[2*i + 1]   *= m_weight / m_dx;
        m_stress[2*i + 1] *= m_weight / m_dx;
    }
    for (int i = 0; i < m_nx + 1; ++i)
    {
        m_mom[2*i + 1] *= m_weight / m_dx;
    }

    // periodic boundary condition for momentum
    m_mom[2*0    + 1] += m_mom[2*m_nx + 1];
    m_mom[2*m_nx + 1] =  m_mom[2*0    + 1];

    // calculate nstress at cell centres and gamma at cell RH faces
    int lh_centre;
    int rh_centre;
    for (int i = 0; i < m_nx; ++i)
    {
        // normalised stress for this cell
        m_nstress[i] = (m_stress[2*i + 0] + m_stress[2*i + 1]) / (m_dens[2*i + 0] + m_dens[2*i + 1]);

        // catch NaN error - this shouldn't have to be used unless there is bad distribution of particles
        if (m_dens[2*i + 0] + m_dens[2*i + 1] == 0)
        {
            m_nstress[i] = 1.0;
            std::cout << "CAUGHT NAN ERROR" << std::endl;
        }
    }
    for (int i = 0; i < m_nx; ++i)
    {
        // calculate necessary indices for RH centre cells for dens and nstress
        if (i < m_nx - 1)
        {
            lh_centre = i;
            rh_centre = i + 1;
        }
        else
        {
            lh_centre = m_nx - 1;
            rh_centre = 0;
        }
        // calculate gamma @ RH face
        // requires: avgmom + mom @ RH face of this cell
        m_gamma[i] = 2.0 * (m_avgmom[i+1] - m_mom[2*(i+1) + 0]) / m_dt;
        // requires: dens + nstress @ this + RH cells

        m_gamma[i] += (0.5 / m_dx) * ((m_dens[2*rh_centre + 0] + m_dens[2*rh_centre + 1]) * m_nstress[rh_centre]
        - (m_dens[2*lh_centre + 0] + m_dens[2*lh_centre + 1]) * m_nstress[lh_centre]);

        m_gamma[i] -= 0.125 * (m_charge / m_mass) * (m_dens[2*lh_centre + 0] + m_dens[2*rh_centre + 0] + m_dens[2*lh_centre + 1] + m_dens[2*rh_centre + 1]) * (m_elec[2*(i + 1) + 0] + m_elec[2*(i + 1) + 1]);

        m_gamma[i] /= 0.5 * (m_dens[2*lh_centre + 1] + m_dens[2*rh_centre + 1]);

        // catch NaN error - this shouldn't have to be used unless there is bad distribution of particles
        if (m_dens[2*lh_centre + 1] + m_dens[2*rh_centre + 1] == 0.0)
        {
            m_gamma[i] = 1.0;
            std::cout << "CAUGHT NAN ERROR" << std::endl;
        }
    }
}


//==============================================================================
/*
push particles with single thread
*/

void Species::push(const double& dt_new)
{
    int p_idx = 83923;

    m_dt = dt_new;

    // set avgmom vector to zeros ahead of the subcycle steps accumulation
    for (int i = 0; i < m_nx + 1; ++i) {
        m_avgmom[i] = 0.0;
    }

    // loop through particles
    for (int p = 0; p < m_np; ++p)
    {
        bool sub_flag  = false;
        int sub_index = -1;
        double subt = 0.0; // subt accumulates the subcycle time that elapses with each subcycle step. at the end of the subcycle, subt = dt

        // transfer current pos, vel and cell to sub-cycle variables
        double subpos  = m_pos[p];
        double subvel  = m_vel[p];
        int subidx = (int) floor(m_pos[p] / m_dx);

        // we track the parent cell idx during the subcycles rather than recalculate it every cycle as otherwise the particle 
        // can get trapped in the cell when we force it to always land at the same face as part of the adaptive push

        // sub-cycle loop
        while (sub_flag == false)
        {
            sub_index++;
            int method_flag = 0;
            double dsubidx = 0;
            double dsubt, dsubpos, dsubvel;

            // get the upper and lower bounds on "dsubpos" -- that is the dsubpos values that push the particle to the lhs or rhs faces of its parent cell
            double dsubpos_lhs_face = (subidx*m_dx) - subpos; // this is the change in position that lands the particle on the lhs face
            double dsubpos_rhs_face = dsubpos_lhs_face + m_dx; // this is the change in position that lands the particle on the rhs face

            // we perform the adaptive push step

            // estimate sub-cycle timestep "subt" with first-order approx, this will propagate the particle approximately "dx" distance
            if (subvel != 0.0) {
                dsubt = m_dx / fabs(subvel);
            } else {
                dsubt = m_dt; 
            }

            // truncate dsubt if it will cause the subcycle elapsed time to exceed "dt"
            if (dsubt + subt > m_dt) {
                dsubt = m_dt - subt;
            }

            // now solve for dsubpos using the Picard method

            // initialise dsubpos with a first order approximation
            dsubpos = dsubt * subvel;

            // ******** solve for dsubpos using Picard method ******** 
            {
                bool pic_flag = false;
                int pic_index = -1;

                double dsubpos_next, dsubvel_next; 

                while (pic_flag == false) {
                    pic_index++;

                    // calculate the electric field at the midstep position (that is subpos + 1/2 * dsubpos)

                    double midstep_pos = subpos + 0.5 * dsubpos; // calculate the midstep position
                  
                    midstep_pos = std::fmod(midstep_pos, m_lx); // we must ensure that midstep_pos lies in [0, lx)
                    if (midstep_pos < 0) {
                        midstep_pos += m_lx;
                    }

                    int midstep_idx = (int) floor(midstep_pos / m_dx); // get the index of the parent cell of the midstep position
                    
                    double midstep_weight_lhs = std::fmax(0, 1.0 - fabs(midstep_pos - midstep_idx*m_dx)/m_dx); // calculate the first-order weights
                    double midstep_weight_rhs = std::fmax(0, 1.0 - fabs(midstep_pos - (midstep_idx + 1)*m_dx)/m_dx);
                    double midstep_elec = 0.5*(m_elec[2*midstep_idx + 0] + m_elec[2*midstep_idx + 1])*midstep_weight_lhs + 0.5*(m_elec[2*(midstep_idx+1) + 0] + m_elec[2*(midstep_idx+1) + 1])*midstep_weight_rhs;

                    // calculate the next iterations of dsubpos and dsubvel using the current iterations of dsubpos and dsubvel
                    // note that the next iteration of dsubpos depends on the next iteration of dsubvel...
                    dsubvel_next = dsubt * (m_charge/m_mass) * midstep_elec;
                    dsubpos_next = dsubt * (subvel + 0.5*dsubvel_next);

                    // if "dsubpos" has converged within the Picard tolerance then end the Picard loop 
                    if (fabs(dsubpos - dsubpos_next) < m_pic_tol) {
                        pic_flag = true;
                    }
                    // if dsubpos hasn't converged in enough Picard iterations then also end the Picard loop
                    if (pic_index > m_max_picard) {
                        pic_flag = true;
                    }

                    // update dsubpos and dsubvel Picard iterations 
                    dsubpos = dsubpos_next;
                    dsubvel = dsubvel_next;
                }
            }

            // if the lhs face of the parent cell has been crossed by the calculated "dsubpos" then we must terminate the particle on this cell face
            // we also alter the parent cell to be the lhs cell
            if (dsubpos <= dsubpos_lhs_face) {
                dsubpos = dsubpos_lhs_face; // force the particle to land on the lhs face
                method_flag = 1; // set the method flag to 1 to indicate that we are forcing the particle to land on the lhs face
                dsubidx = -1; // decrement the parent cell idx to show that the particle is now in the lhs cell
            } else if (dsubpos >= dsubpos_rhs_face) {
                dsubpos = dsubpos_rhs_face; // force the particle to land on the rhs face
                method_flag = 1; // set the method flag to 1 as we force the particle to land on rhs face
                dsubidx = 1; // increment the parent cell idx as the particle is now in the lhs cell
            }
            // we must be sure to reset method_flag and dsubidx to be zero before the next subcycle!

            // if the lhs or rhs cell face has been crossed then we must recalculate the "dsubt" subcycle time step that elapses for this "dsubpos"
            if (method_flag == 1) {
                // calculate the electric field at the midstep for the "dsubpos" that we have now calculated
                int midstep_idx = (int) floor((subpos + 0.5*dsubpos) / m_dx);
                double midstep_weight_lhs = 1.0 - (subpos+0.5*dsubpos - midstep_idx*m_dx)/m_dx;
                double midstep_weight_rhs = 1.0 - midstep_weight_lhs;
                double midstep_elec = 0.5*(m_elec[2*midstep_idx + 0] + m_elec[2*midstep_idx + 1])*midstep_weight_lhs + 0.5*(m_elec[2*(midstep_idx+1) + 0] + m_elec[2*(midstep_idx+1) + 1])*midstep_weight_rhs;

                // ******** solve for "dsubt" using the Picard method ******** 
                // we already have an estimate for "dsubt", so the method can run without needing to initialise "dsubt"
                {
                    bool pic_flag = false;
                    int pic_idx = -1;

                    double dsubt_next;

                    while (pic_flag == false) {
                        pic_idx++;

                        // calculate the next increment of "dsubt" by combining both equations for "dsubpos" and "dsubvel"
                        dsubt_next = dsubpos / (subvel + 0.5*dsubt*(m_charge/m_mass)*midstep_elec); // there is a potential for NAN error here, perhaps guard against this

                        // stop the estimate from being negative
                        if (dsubt_next < 0.0) {
                            dsubt_next = 0.0;
                        }
                        
                        // if "dsubt" has converged within the Picard tolerance then end the Picard loop 
                        if (fabs(dsubt - dsubt_next) < m_pic_tol) {
                            pic_flag = true;
                        }
                        // if dsubt hasn't converged in enough Picard iterations then also end the Picard loop
                        if (pic_idx > m_max_picard) {
                            pic_flag = true;
                        }

                        // update the Picard iteration
                        dsubt = dsubt_next;
                    }
                }

                // calculate "dsubvel" using the "dsubpos" and the result we obtained for "dsubt"
                dsubvel = dsubt * (m_charge/m_mass) * midstep_elec;
            }

            // at this stage we now have "dsubpos", "dsubvel" and "dsubt" values calculated
            // we check if this value of "dsubt" takes us up to the require elapsed time for the subcycles
            if ((subt + dsubt >= m_dt) || (sub_index > 10000))
                sub_flag = true;

            // accumulate avgmom for this subcycle step
            // the algorithm has ensured that the particle has not travelled further than the lhs and rhs faces of its parent cell
            // this means that for this subcycle step we accumulate momentum onto the lhs and rhs faces of the parent cell, given by "subidx"
            // first calculate the first-order weights for "dsubpos" to interpolate the electric field, we already know that the midstep_idx = subidx
            double midstep_weight_lhs = 1.0 - ((subpos + 0.5*dsubpos) - subidx*m_dx)/m_dx;
            double midstep_weight_rhs = 1.0 - midstep_weight_lhs;
            m_avgmom[subidx]     += (1.0 / (m_dx * m_dt)) * m_weight * (subvel + 0.5*dsubvel) * midstep_weight_lhs * dsubt;
            m_avgmom[subidx + 1] += (1.0 / (m_dx * m_dt)) * m_weight * (subvel + 0.5*dsubvel) * midstep_weight_rhs * dsubt;

            if (p == p_idx) {
                std::cout << "subcycle index = " << sub_index << " -- dsubpos = " << dsubpos << ", dsubt = " << dsubt << std::endl;
                std::cout << "dt = " << m_dt << ", subt + dsubt = " << subt + dsubt << std::endl;
            }

            // now we calculate the new "subpos", "subidx" values and apply boundary conditions
            double subpos_new = subpos + dsubpos;
            int subidx_new = subidx + dsubidx;
            if (subidx_new < 0) {
                subidx_new += m_nx;
                subpos_new += m_lx;
            } else if (subidx_new > (m_nx - 1)) {
                subidx_new -= m_nx;
                subpos_new -= m_lx;
            }

            subt += dsubt;
            subvel += dsubvel;
            subpos = subpos_new;
            subidx = subidx_new;
        }

        // transfer sub-cycle new values to global new values at end of particle evolution
        m_pos_new[p] = subpos;
        m_vel_new[p] = subvel;
    }

    // apply periodic boundary conditions to average momentum 
    m_avgmom[0] += m_avgmom[m_nx];
    m_avgmom[m_nx] = m_avgmom[0];
}

//==============================================================================
//==============================================================================


/*
step HO values - pos, vel, cell, dens and stress -> shift future values to current position and then extrapolate new current values to new future position
*/
void Species::step()
{
    for (int i = 0; i < m_np; ++i) {
        m_pos[i] = m_pos_new[i];
        m_vel[i] = m_vel_new[i];
    }

    // project future values to new current values for density and stress
    for (int i = 0; i < m_nx; ++i)
    {
        m_dens[2*i + 0]   = m_dens[2*i + 1];
        m_stress[2*i + 0] = m_stress[2*i + 1];
    }
    for (int i = 0; i < m_nx + 1; ++i)
    {
        m_mom[2*i + 0] = m_mom[2*i + 1];
    }
}

// getters
double* Species::get_dens_ptr()
{
    return m_dens;
}

double* Species::get_mom_ptr()
{
    return m_mom;
}

double* Species::get_avgmom_ptr()
{
    return m_avgmom;
}

double* Species::get_nstress_ptr()
{
    return m_nstress;
}

double* Species::get_gamma_ptr()
{
    return m_gamma;
}

void Species::clear_files()
{
    std::vector<std::string> files = {
        "continuity", "dens", "mom", "avgmom", "stress", "gamma", "pos", "vel"
    };

    for (const auto &file : files) {
        std::ofstream outfile("./output/ho_" + m_name + "_" + file + ".csv", std::ofstream::out | std::ofstream::trunc);
        outfile.close();
    }
}

void print_array_csv(std::string filepath, double *data, int N, int stride = 1, int offset = 0) {
    // open the file
    std::ofstream file(filepath, std::ofstream::out | std::ofstream::app);
    
    // print N values from the array
    for (int i = 0; i < N; ++i) {
        file << data[i*stride + offset];
        if (i < N - 1) {
            // only print "," if not the last value on the line
            file << ",";
        }
    }
    file << "\n";

    // close the file
    file.close();
}

void Species::output_data() {
    // do not calculate continuity for initial timestep as there is no valid avgmom
    // analysis data
    /*for (int i = 0; i < nx; ++i)
    {
        continuity_res_single[i] = 0.0;
    }

    for (int i = 0; i < nx; ++i)
    {
        continuity_res_single[i] = (dens_single[2*i + 1] - dens_single[2*i + 0])/dt;
        continuity_res_single[i] += (avgmom_single[i + 1] - avgmom_single[i])/dx;
    }*/

    // continuity
    //print_array_csv("./output/ho_" + name + "_continuity.csv", continuity_res_single, nx, 1, 0);

    // density
    print_array_csv("./output/ho_" + m_name + "_dens.csv", m_dens, m_nx, 2, 1);

    // momentum
    print_array_csv("./output/ho_" + m_name + "_mom.csv", m_mom, m_nx + 1, 2, 1);

    // avg momentum
    print_array_csv("./output/ho_" + m_name + "_avgmom.csv", m_avgmom, m_nx + 1, 1, 0);

    // stress
    print_array_csv("./output/ho_" + m_name + "_stress.csv", m_stress, m_nx, 2, 1);

    // gamma
    print_array_csv("./output/ho_" + m_name + "_gamma.csv", m_gamma, m_nx, 1, 0);
}