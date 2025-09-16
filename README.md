# Vlasov-Ampère solvers

![Build](https://img.shields.io/github/actions/workflow/status/username/repo/ci.yml?branch=main)
![License](https://img.shields.io/github/license/username/repo)
![Contributors](https://img.shields.io/github/contributors/username/repo)
#### Authorship
The solvers implemented in this project are improved versions of those implemented in the 2022 Final Year Project for MEng Aeronautical Engineering, Imperial College by William Denny. The background details provided here have been have been lifted from the project report. 

## Table of Contents
- [Background](#background)
- [Usage](#usage)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
- [Requirements](#requirements)

<!-- <!-- mkdir build && cd build -->

## Background
The Vlasov-Amp&egrave;re system of equations can be solved to simulate the behaviour of plasma. This project is concerned with implementing System 2.4 of the Project NEPTUNE specification: solving the Vlasov-Amp&egrave;re system of equations using the kinetic-enslaved particle-in-cell (PIC) implicit moment method (IMM). The Vlasov-Amp&egrave;re equations characterise the velocity distribution of plasma charge carriers, which are solved using the PIC method by simulating charge carriers as 'super particles'.

### Theory

We provide some technical details about electromagnetism that are relevant in the construction of the kinetic-enslaved PIC IMM model. In a plasma of charged particles coupled with an electromagnetic field, a single particle will experience a Lorentz force from the electrostatic and magnetic components of field, given by

$$ \mathbf{F} = q_p \mathbf{E} + q_p(\mathbf{v_p} \times \mathbf{B}), \qquad (1) $$

where $\mathbf{F}$ is the force vector, $\mathbf{E}$ is the electric field vector, $\mathbf{B}$ is the magnetic field vector, $\mathbf{v_p}$ is the velocity vector $q_p$ is the particle charge. One consequence of (1) is that the magnetic field will not exert a force on a particle with velocity
 parallel to its field lines. In tokamak fusion reactors, a magnetic field is used to direct plasma flow
 around a circular path. In equilibrium, the plasma flow is parallel to the field lines, meaning the force
 exerted by the magnetic field can be neglected [2].

Maxwell's equations describe the dynamics of electromagnetic fields using the moment quantities of charged particles, that is, charge density $\rho$ and current vector $\mathbf{j}$. These are given by

$$ \nabla \cdot \mathbf{E} = \frac{\rho}{\epsilon_0} \quad (2) \qquad \nabla \times \mathbf{E} = -\mu_0 \frac{\partial \mathbf{B}}{\partial t} \quad (3) $$
$$ \nabla \cdot \mathbf{B} = 0 \quad (4) \qquad \nabla \times \mathbf{B} = -\mu_0 \mathbf{j} + \mu_0 \epsilon_0 \frac{\partial \mathbf{E}}{\partial t} \quad (5) $$

The electric field $\mathbf{E}$ can be expressed as the gradient vector field of the electric potential $\phi$, in equation (6). (2) can then be written in terms of $\phi$ to get the Poisson equation (7):

$$\mathbf{E} = -\nabla \phi \quad (6) \qquad \nabla^2 \phi = -\frac{\rho}{\epsilon_0} \quad (7)$$

(4) is known as the Amp&egrave;re equation. By neglecting magnetic field, the Amp&egrave;re the equation can be expressed by (8).

$$ \epsilon_0 \frac{\partial \mathbf{E}}{\partial t} + \mathbf{j} - \langle \mathbf{j} \rangle \quad (8) \qquad \mathbf{j} = q_\alpha \mathbf{nu} \quad (9)$$

$\mathbf{nu}$ is the charge-normalised current, also referred to in [1] as the number-momentum or momentum. Moreover $n$ is the charge-normalised density, that is, $\rho = q n$, and is referred to in [1] as number density or density.

The velocity distrbution function of a plasma $f(\mathbf{x},\mathbf{v},t)$ is the probability distribution function over $\mathbb{R}^3$ in the velocity space. More precisely, $f$ quantifies the probability that a particle at location $\mathbf{x}$ at time $t$ has velocity $\mathbf{v}$. Most gases usually have a Maxwell-Boltzmann velocity distribution that is invariant with time. However, plasmas coupled to electromagnetic fields have a velocity distribution that evolves with time according to the Vlasov equation. In the case that we neglect the magnetic field, the Vlasov equation in 1D is:

$$\frac{\partial f}{\partial t} + v \frac{\partial f}{\partial x} + \frac{q}{m} E \frac{\partial f}{\partial v}, \qquad (10)$$

where $m$ is the particle mass. Taking the zeroth and first moments of the Vlasov equation returns the continuity and momentum conservation equations:

$$ \frac{\partial n}{\partial t} + \frac{\partial (nu)}{\partial x} = 0 \quad (11) \qquad 
\frac{\partial (nu)}{\partial t} + \frac{\partial S}{\partial x} - \frac{q}{m} n E \quad (12)

$$

(11) and (12) are then combined with (8) to get the Vlasov-Amp&egrave;re system of equations as given in [1].

The PIC model represents the velocity distribution as a collection of super particles. In [1], there is an lower-order (LO) system of density, momentum and electric field discretised over a mesh of cells. The higher-order (HO) system is the collection of super particles. The LO system is solved implicitly to obtain the electric field and estimates of the density and momentum. The HO system is solved by integrating the Vlasov equation in the Lagrangian frame of reference for each super particle $p$:

$$\frac{\partial v_p}{\partial t} = \frac{q}{m} E_p \quad (13) \qquad \frac{\partial x_p}{\partial t} = v_p \quad (14) $$


 The particle moments are accumulated to solve for the electromagnetic field which is discretised over a mesh of cells.

 Kinetic-enslavement is the process of coupling the LO and HO systems together. [1] provides details on how this is done. It mostly involves include accumulating the momentum and density from HO system and injecting this into the LO system equations. 

## Usage
To begin simply clone the repository
```bash
git clone https://github.com/greenvale/vlasov-ampere.git
cd vlasov-ampere
```
### Explicit solver
The explicit solver is implemented by using a ``Species`` class to encapsulate the particle species and solve for their density, momentum, position and velocity. The electric field is then solved seperately and is propagated to each species during each time step.

```cpp
#include "library.h"
int main() {
    simulate();
}
```

## Requirements

This project requires an NVidia GPU with compute capability &ge; 6.0. CUDA Toolkit &ge; 12.0 must also be installed to run the ``nvcc`` compiler.

## Demo
![Simulation gif](docs/demo.gif)

### References

[1] Smith, J., et al., *A Fast CUDA Particle Simulator*, Journal of Computational Physics, 2020. [Link](https://doi.org/10.1000/exampledoi)  
[2] Doe, A., *GPU-Accelerated Simulations*, Conference on HPC, 2019. [Link](https://doi.org/10.1001/exampledoi)
