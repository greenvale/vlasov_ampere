#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <memory>
#include <utility>
#include <vector>
#include <random>
#include <fstream>
#include <cuda_runtime.h>
#include <numeric>
#include <cmath>
#include <functional>
#include <cstring>

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

namespace constants
{
    const double e0 = 1.0;//8.85418 * pow(10.0, -12);
    const double Kb = 1.0;//1.3806 * pow(10.0, -23);
}

//**********************************************************************************************************************************
// Kernels


__global__ void n_kernel(double *x, double *n, int Np, int Ncells, double length) {
    // create a local buffer for the block
    extern __shared__ double n_local[];

    // initialise the local number density
    if (threadIdx.x < Ncells + 1) {
        n_local[threadIdx.x] = 0.0;
    }
    
    __syncthreads();

    double dx = length / Ncells;

    // get the global index for this thread among all threads in all blocks
    int idx = threadIdx.x + blockDim.x * blockIdx.x;

    for ( ; idx < Np; idx += blockDim.x * gridDim.x) {
        // get the index of the lhs cell face
        int lhs_idx = int(x[idx] / dx);

        //if (lhs_idx < 0 || lhs_idx >= Ncells) {
        //    printf("x=%f idx=%d lhs_idx=%d\n", x[idx], idx, lhs_idx);
        //}

        // assert that the lhs_idx is in range, i.e. x lies in [0, length)
        // this shouldn't be necessary to check but is important to catch if it does go wrong
        assert(lhs_idx >= 0 && lhs_idx < Ncells);

        // now get the shape function values -- the weights that the particle exert
        // lhs_weight is the weight exerted by the particle at the lhs cell face, and vice versa for rhs_weight
        double lhs_weight = 1.0 - (x[idx] - lhs_idx*dx) / dx;
        double rhs_weight = 1.0 - lhs_weight;

        // accumulate the weights on the local buffers
        atomicAdd(&n_local[lhs_idx], lhs_weight);
        atomicAdd(&n_local[lhs_idx + 1], rhs_weight);
    }

    __syncthreads();

    // assign (Ncells + 1) many threads in this block to accumulate the local buffer for this block
    // into the global buffer for the density
    if (threadIdx.x < Ncells + 1) {
        atomicAdd(&n[threadIdx.x], n_local[threadIdx.x]);
    }
}

__global__ void nu_kernel(double *x, double *v, double *nu, int Np, int Ncells, double length) {
    // create a local buffer for the block
    extern __shared__ double nu_local[];

    // initialise the local momentum buffer using (Ncells + 1) many threads in this block
    if (threadIdx.x < Ncells + 1) {
        nu_local[threadIdx.x] = 0.0;
    }
    
    __syncthreads();

    double dx = length / Ncells;

    int idx = threadIdx.x + blockDim.x * blockIdx.x;

    for ( ; idx < Np; idx += blockDim.x * gridDim.x) {
        int lhs_idx = int(x[idx] / dx);

        //if (lhs_idx < 0 || lhs_idx >= Ncells) {
        //    printf("x=%f idx=%d lhs_idx=%d\n", x[idx], idx, lhs_idx);
        //}

        assert(lhs_idx >= 0 && lhs_idx < Ncells);

        double lhs_weight = 1.0 - (x[idx] - lhs_idx*dx) / dx;
        double rhs_weight = 1.0 - lhs_weight;

        atomicAdd(&nu_local[lhs_idx], lhs_weight * v[idx]);
        atomicAdd(&nu_local[lhs_idx + 1], rhs_weight * v[idx]);
    }

    __syncthreads();

    if (threadIdx.x < Ncells + 1) {
        atomicAdd(&nu[threadIdx.x], nu_local[threadIdx.x]);
    }
}

__global__ void xv_kernel(double *x, double *v, double *E, int Np, int Ncells, double length, double dt, double charge, double mass) {

    int idx = threadIdx.x + blockDim.x * blockIdx.x;

    double dx = length / Ncells;

    for ( ; idx < Np; idx += blockDim.x * gridDim.x) {
        int lhs_idx = int(x[idx] / dx);

        //if (lhs_idx < 0 || lhs_idx >= Ncells) {
        //    printf("x=%f idx=%d lhs_idx=%d\n", x[idx], idx, lhs_idx);
        //}

        assert(lhs_idx >= 0 && lhs_idx < Ncells);

        // calculate the lhs and rhs weights for this particle w.r.t. the cell faces
        double lhs_weight = 1.0 - (x[idx] - lhs_idx*dx) / dx;
        double rhs_weight = 1.0 - lhs_weight;

        // calculate the acceleration experienced by the particle by interpolating the electric field using the weights
        double accel = (charge / mass) * (E[lhs_idx]*lhs_weight + E[lhs_idx + 1]*rhs_weight);

        // calculate the change in position and velocity
        // the position change is second order explicit, the velocity change is first order explicit
        double dx = v[idx]*dt + 0.5*accel*dt*dt;
        double dv = accel * dt;

        // update the position and velocity
        x[idx] += dx;
        v[idx] += dv;

        // apply periodic boundary conditions to position and velocity
        x[idx] = fmod(x[idx], length);
        if (x[idx] < 0) {
            x[idx] += length;
        }
    }
}


//**********************************************************************************************************************************
// Helper functions

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


//**********************************************************************************************************************************
// Species class

class Species {
    public:
        double *x_d;
        double *v_d;
        double *n_h, *n_d;
        double *nu_h, *nu_d;
        double weight;
        double charge;
        double mass;
        int Np;
        
        // domain information
        int Ncells;
        double length;

        // counters that keep track of the number of times the x, v, n or nu values have been upated
        int x_counter = 0;
        int v_counter = 0;
        int n_counter = 0;
        int nu_counter = 0;

        bool initialised = false;

    Species() {

    }

    //void initialise(int Np, int Ncells, int Nvcells, double length, double v_range, double charge, double target_n, std::vector<double> pdf, std::mt19937 &gen) {
    void initialise(int Np, int Ncells, int Nvcells, double length, double v_range, double charge, double mass, double n_avg_init, 
        std::function<double(double,double,double,double)> unscaled_pdf, std::mt19937 &generator) {
        if (initialised == false) {
            // copy to member variables
            this->charge = charge;
            this->length = length;
            this->mass = mass;
            this->Np = Np;
            this->Ncells = Ncells;
            
            // initialise the host memory
            n_h = new double[Ncells + 1];
            nu_h = new double[Ncells + 1];

            // initialise the device memory
            CUDA_CHECK(cudaMalloc((void**) &x_d, Np * sizeof(double)));
            CUDA_CHECK(cudaMalloc((void**) &v_d, Np * sizeof(double)));
            CUDA_CHECK(cudaMalloc((void**) &n_d, (Ncells + 1) * sizeof(double)));
            CUDA_CHECK(cudaMalloc((void**) &nu_d, (Ncells + 1) * sizeof(double)));

            // initialise the particle position and velocity using the pdf
            // the particles will be initialised on the host and then the pos/vel data
            // will be copied across to the device memory
            double dx = length / Ncells;
            double dv = 2 * v_range / Nvcells;

            double *x_init, *v_init;

            x_init = new double[Np];
            v_init = new double[Np];

            // create a uniform distribution object which be sampled from to distribute particles
            std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

            // get the pdf grid -- that is a grid of weightings for each x-v cell that sum to 1 and indicate
            // how many particles will be created in each cell
            std::vector<double> pdf_grid = get_pdf_grid(Ncells, Nvcells, length, v_range, unscaled_pdf);
            
            int ptcl_idx = 0;
            for (int i = 0; i < Ncells; ++i) {
                for (int j = 0; j < Nvcells; ++j) {
                    // calculate the target number of particles for this x-v cell
                    int target_np = int(Np * pdf_grid[i*Nvcells + j]);
                    for (int k = 0; k < target_np; ++k) {
                        // particles are uniformly distributed in this x-v cell
                        double x = (i + uniform_distribution(generator)) * dx;
                        double v = -v_range + (j + uniform_distribution(generator)) * dv;

                        x_init[ptcl_idx] = x;
                        v_init[ptcl_idx] = v;
                        ptcl_idx++;
                    }
                }
            }

            // randomly place the remaining particles that need to be initialised
            for ( ; ptcl_idx < Np; ++ptcl_idx) {
                // generate a random x,v value
                double x = uniform_distribution(generator) * length;
                double v = -v_range + uniform_distribution(generator) * 2 * v_range;
                x_init[ptcl_idx] = x;
                v_init[ptcl_idx] = v;
            }

            // now copy the initial particle positions and velocities from host to device
            CUDA_CHECK(cudaMemcpy(x_d, x_init, Np * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(v_d, v_init, Np * sizeof(double), cudaMemcpyHostToDevice));

            // delete the initial position and velocity data from the host memory
            delete[] x_init;
            delete[] v_init;

            // now we calculate the "weight" of each particle. that is the weighting it gets
            // when calculating the density and momentum. We first assume that the weight is 1
            // for each particle, calculate the density, then calculate the weight such that
            // the average density is equal to the target average density.

            // set the density to zeros before accumulating
            CUDA_CHECK(cudaMemset(n_d, 0, (Ncells + 1) * sizeof(double)));

            // accumulate the density. The density kernel requires that the position of each particle
            // lies in [0, length). There are no checks in the density kernel as they should not be required.
            int blockSize = 265;
            int numBlocks = 32;
            n_kernel<<<numBlocks, blockSize, (Ncells + 1) * sizeof(double)>>>(x_d, n_d, Np, Ncells, length);

            // wait for the kernel to finish
            CUDA_CHECK(cudaDeviceSynchronize());
            
            // copy density from device to host
            CUDA_CHECK(cudaMemcpy(n_h, n_d, (Ncells + 1) * sizeof(double), cudaMemcpyDeviceToHost));

            // apply boundary conditions to density
            n_h[0] += n_h[Ncells];
            n_h[Ncells] = n_h[0];

            // calculate the average density when weight = 1.0
            // you only accumulate across Ncells as the far rhs cell face is the far lhs cell face
            double n_avg = std::accumulate(n_h, n_h + Ncells, 0.0);
            n_avg /= Ncells;
            this->weight = n_avg_init / n_avg;

            std::cout << "Weight = " << weight << std::endl;

            // mark that this species has been initialised meaning this function cannot be run again
            initialised = true;
        }
    }

    ~Species() {
        if (this->initialised) {
            cudaFree(x_d);
            cudaFree(v_d);
            cudaFree(n_d);
            cudaFree(nu_d);
            delete[] n_h;
            delete[] nu_h;
        }
    }

    void update_n(int blockSize, int numBlocks) {
        // set the density to zeros before accumulating
        CUDA_CHECK(cudaMemset(n_d, 0, (Ncells + 1) * sizeof(double)));

        // accumulate the density
        n_kernel<<<numBlocks, blockSize, (Ncells + 1) * sizeof(double)>>>(x_d, n_d, Np, Ncells, length);

        CUDA_CHECK(cudaDeviceSynchronize());

        // copy density from device to host
        CUDA_CHECK(cudaMemcpy(n_h, n_d, (Ncells + 1) * sizeof(double), cudaMemcpyDeviceToHost));

        // scale the density by the particle weight, in the current calculation each particle has weight = 1
        for (int i = 0; i < Ncells + 1; ++i) {
            n_h[i] *= this->weight;
        }

        // apply bounary conditions 
        n_h[0] += n_h[Ncells];
        n_h[Ncells] = n_h[0];

        // increment the n counter
        n_counter++;
    }

    void update_nu(int blockSize, int numBlocks) {
        // set the momentum to zeros before accumulating
        CUDA_CHECK(cudaMemset(nu_d, 0, (Ncells + 1) * sizeof(double)));

        // accumulate the momentum
        nu_kernel<<<numBlocks, blockSize, (Ncells + 1) * sizeof(double)>>>(x_d, v_d, nu_d, Np, Ncells, length);

        // wait for the kernel to finish
        CUDA_CHECK(cudaDeviceSynchronize());

        // copy the calculated momentum from the device to host
        CUDA_CHECK(cudaMemcpy(nu_h, nu_d, (Ncells + 1) * sizeof(double), cudaMemcpyDeviceToHost));

        // scale the momentum by the particle weight, in the current calculation each particle has weight = 1
        for (int i = 0; i < Ncells + 1; ++i) {
            nu_h[i] *= this->weight;
        }

        // apply the boundary conditions
        nu_h[0] += nu_h[Ncells];
        nu_h[Ncells] = nu_h[0];

        // increment the nu counter
        nu_counter++;
    }

    void update_xv(int blockSize, int numBlocks, double *E_h, double dt) {
        // E_h is a pointer to the electric field buffer in the host memory
        // copy E_h into the device memory
        double *E_d;
        CUDA_CHECK(cudaMalloc((void**) &E_d, (Ncells + 1) * sizeof(double)));
        CUDA_CHECK(cudaMemcpy(E_d, E_h, (Ncells + 1) * sizeof(double), cudaMemcpyHostToDevice));

        // push the particles
        xv_kernel<<<numBlocks, blockSize>>>(x_d, v_d, E_d, Np, Ncells, length, dt, this->charge, this->mass);

        // wait for the kernel to finish
        CUDA_CHECK(cudaDeviceSynchronize());

        // now free the temporary electric field data on the device
        cudaFree(E_d);

        x_counter++;
        v_counter++;
    }

};


double unscaled_pdf_landau(double length, double v_range, double x, double v) {
    double fx = 1 + (0.2 * sin(2.0 * M_PI * x / length));
    double fv = std::exp(-v*v / (2 * (v_range*v_range / 1.96*1.96)));
    return fx * fv;
}

int main(int argc, char *argv[])
{
    int device;
    cudaGetDevice(&device);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device);
    std::cout << "Device name: " << prop.name << "\n";
    std::cout << "Compute capability: " 
              << prop.major << "." << prop.minor << "\n";
              
    // open data file
    std::ofstream n_outfile("n.csv", std::ios::trunc);
    std::ofstream nu_outfile("nu.csv", std::ios::trunc);
    std::ofstream E_outfile("E.csv", std::ios::trunc);

    // random generator
    std::random_device rd;
    std::mt19937 generator(rd());
    
    // parameters
    int Ncells = 16;
    int Nvcells = 8;
    int Np = 2 << 20;
    double dt = 0.0005;
    double length = 1.0;
    double v_range = 8.0;
    double target_n = 1.0;

    // calculated parameters
    double dx = length / Ncells;

    // declare pointers for electric field
    double *E_h;//, *E_new_h;

    // device parameters
    int blockSize = 265;
    int numBlocks = 32;

    // allocate host buffer for electric field
    E_h = new double[Ncells + 1];
    //E_new_h = new double[Ncells + 1];

    // create electrons species
    Species *p_electrons = new Species;

    // initialise the species
    p_electrons->initialise(Np, Ncells, Nvcells, length, v_range, -1.0, 1.0, target_n, unscaled_pdf_landau, generator);
    
    // Now initialise the electric field
    // calculate the density of each species, then calculate the overall charge density
    // then calculate the electric field by solving the poisson equation.

    // update the density for the species
    p_electrons->update_n(blockSize, numBlocks);

    // calculate the charge density - this is the n-density multiplied by the species charge
    // then summed over all species to get the overall charge density
    double *charge_n_init = new double[Ncells + 1];

    for (int i = 0; i < Ncells + 1; ++i) {
        charge_n_init[i] = p_electrons->n_h[i] * p_electrons->charge;
    }
    
    // assume the presence of ions with charge 1
    for (int i = 0; i < Ncells + 1; ++i) {
        charge_n_init[i] += 1.0;
    }

    // calculate the initial electric field using the charge density
    // dE/dx = charge_dens / e0
    // => E_{i+1} = E_{i} + dx * charge_dens_{i} / e0
    E_h[0] = 0.0;
    for (int i = 0; i < Ncells; ++i) {
        E_h[i + 1] = E_h[i] + charge_n_init[i] * dx / constants::e0; 
    }
    //E_h[Ncells] = E_h[Ncells - 1] + E_h[1] - E_h[0];
    E_h[Ncells] = 0.0;

    // Now calibrate the electric field -- calculate the average value and then
    // subtract from E to ensure that the integral is zero overall.
    double E_avg = std::accumulate(E_h, E_h + Ncells, 0.0);
    E_avg /= Ncells;
    for (int i = 0; i < Ncells + 1; ++i) {
        E_h[i] -= E_avg;
    }

    for (int k = 0; k < 5000; ++k) {
        
        //****************************************************************************************
        // print values
        if (k % 10 == 0) {
            // save density data
            if (n_outfile.is_open()) {
                for (int i = 0; i < Ncells + 1; ++i) {
                    n_outfile << p_electrons->n_h[i];
                    if (i < Ncells) {
                        n_outfile << ",";
                    }
                }
                n_outfile << "\n";
            }

            // save electric field data
            if (E_outfile.is_open()) {
                for (int i = 0; i < Ncells + 1; ++i) {
                    E_outfile << E_h[i];
                    if (i < Ncells) {
                        E_outfile << ",";
                    }
                }
                E_outfile << "\n";
            }

            // save momentum data
            if (nu_outfile.is_open()) {
                for (int i = 0; i < Ncells + 1; ++i) {
                    nu_outfile << p_electrons->nu_h[i];
                    if (i < Ncells) {
                        nu_outfile << ",";
                    }
                }
                nu_outfile << "\n";
            }
        }
        //****************************************************************************************
        
        // calculate the density for each species
        p_electrons->update_n(blockSize, numBlocks);

        // calculate the momentum for each species 
        p_electrons->update_nu(blockSize, numBlocks);

        // push particles through the current electric field to get new position, velocity data
        p_electrons->update_xv(blockSize, numBlocks, E_h, dt);

        // calculate the temporal change in electric field using the current momentum from the current velocity data
        // set the new electric field to zero
        //std::fill(E_new_h, E_new_h + Ncells + 1, 0.0);

        // get the average charged momentum across all species 
        double nu_avg_electrons = p_electrons->charge * std::accumulate(p_electrons->nu_h, p_electrons->nu_h + Ncells + 1, 0.0);
        nu_avg_electrons /= Ncells;

        // calculate the new electric field using dE/dt = q/e0 * nu - \bar{nu}
        for (int i = 0; i < Ncells + 1; ++i) {
            E_h[i] += (dt / constants::e0) * (p_electrons->nu_h[i] - nu_avg_electrons);
        }

        // update the electric field to be the new electric field
        //std::memcpy(E_h, E_new_h, (Ncells + 1) * sizeof(double));

    }
    
}