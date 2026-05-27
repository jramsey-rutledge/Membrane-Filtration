# Membrane-Filtration Model
This project implements a two-dimensional finite volume solver for the steady advection–diffusion equation to model concentration accumulation in a membrane filtration system. The solver investigates transport behavior driven by bulk flow advection, molecular diffusion, and membrane permeation effects.

The code supports both:
 - uniform Cartesian grids, and
 - nonuniform stretched meshes for enhanced near-wall resolution.

The governing equations are discretized using the finite volume method (FVM) and assembled into a sparse linear system solved using the Eigen C++ library.

**Author:**
- Jenna Ramsey-Rutledge

## Physical Background
Membrane filtration systems are widely used in water purification, desalination, and biomedical separation processes. During filtration, solute particles are transported toward the membrane surface by advection while diffusion acts to redistribute concentration gradients.

As permeate is extracted through the membrane, solute accumulation can occur near the membrane wall, leading to concentration polarization. This localized increase in concentration reduces filtration efficiency and can contribute to membrane fouling.

This project models the steady transport of a passive scalar concentration field within a channel subject to:
 - streamwise advection,
 - transverse membrane permeation,
 - and molecular diffusion.

The solver is intended as a simplified CFD framework for studying concentration accumulation and near-wall transport phenomena in membrane systems.

## Physics and Mathematical Model

### Governing Equations
The concentration field is governed by the steady two-dimensional advection–diffusion equation:

$$ u \frac{\partial C}{\partial x} + v \frac{\partial C}{\partial y} = D \left( \frac{\partial^2 C}{\partial x^2} + \frac{\partial^2 C}{\partial y^2} \right)$$

where $C$ is the scalar concentration, $u$ and $v$ are the velocity components, and $D$ is the diffusion coefficient. The equation is solved in conservative finite-volume form using face fluxes and diffusive transport coefficients.

### Boundary Conditions
The membrane walls are defined on the North and South boundaries, while the West boundary is the inlet and East boundary is the outlet. The membrane walls are modeled using Robin boundary conditions. The North wall is defined with 

$$ -V_m C + D \frac{\partial C}{\partial y} = 0 $$

and the South wall is defined with 

$$ V_m C - D \frac{\partial C}{\partial y} = 0 $$

where $V_m$ is the membrane permeation velocity. These conditions model the balance between membrane suction and diffusive transport at the membrane surface.

The West boundary utilizes a fixed inlet concentration perscribed as $C=1$. This is implemented as a Dirichlet boundary contion. The East boundary utilizes a zero-gradient outflow condition applied as $\frac{\partial C}{\partial x} = 0$. The outlet boundary is implemented using a second-order finite difference extrapolation stencil. 

## Numerical Methods 
### Spatial Discretization

The governing equations are discretized using the finite volume method on a structured Cartesian mesh. Each control volume contains east/west diffusive fluxes, north/south diffusive fluxes, and convective transport through cell faces.A five-point stencil is used for the interior control volumes.

### Convective Fluxes

Convective transport terms are discretized using central differencing. Face fluxes are computed separately for east/west faces and north/south faces. The discretized convective contribution is assembled directly into the sparse coefficient matrix.

### Diffusion Terms

Diffusive transport coefficients are computed using face areas and local cell spacing:

\begin{align}
$$a_E & = \frac{D A_E}{\Delta x_E}$$\\
$$a_N & = \frac{D A_N}{\Delta y_N}$$
\end{align}

with analogous expressions for west and south faces.

### Nonuniform Mesh Clustering
Although the solver supports both uniform mesh spacing and nonuniform stretched spacing, the nonuniform mesh is best suited for this setting because it improves near-wall resolution and better captures the steep concentration gradients associated with concentration polarization. The nonuniform mesh generator uses hyperbolic tangent stretching functions to cluster grid points near the membrane walls. The wall-normal face locations are generated using:

$$ y_f = \frac{L_y}{2} \left( 1 + \frac{\tanh(\beta(2s-1))}{\tanh(\beta)}\right) $$ 

where $s \in \[0,1\]$ is the normalized computational coordinate and $\beta$ controls clustering intensity.

### Linear System Assembly and Solver

The discretized equations are assembled into a sparse linear system:

$$A\phi=b$$

where $A$ is the sparse coefficient matrix, $\phi$ is the concentration vector, and $b$ is the right-hand side vector. Sparse matrix coefficients are stored using triplet assembly before constructing the Eigen sparse matrix object.

The resulting linear system is solved using the Eigen library sparse LU decomposition solver. The solution is then mapped back into the scalar field storage structure for post-processing and visualization.


## Quick Start

### Prerequisites

- C++17-compatible compiler
- GNU Make
- Eigen 5.0.0+

### Building and Running

```bash
# Clone the repository
git clone <repo-url>
cd Membrane-Filtration

# Change eigen path in Makefile

# Compile the project
make build

# Run the program
./build/solver.exe

# Remove compiled artifacts
make clean

# Build and run in one step
make
```

### Adjusting the Simulation
Edit the .inp files in 'input/' to change the spatial domain, temporal domain, fan parameters, or physical properties before building

## Project Structure

```
Membrane-Filtration/
├── src/                      # Source and Header Files
│   ├── main.cpp
│   ├── matrixType.{cpp,h}           # Linear System
│   ├── meshType.{h}                 # Spatial Mesh Parent Class
│   ├── nonUniformGrid.{cpp,h}       # Non-Uniform Mesh
│   ├── uniform Grid.{cpp,h}         # Uniform Mesh
│   ├── simTime.{cpp,h}              # Time Stepping  
│   ├── volScalarField.{cpp,h}       # Temperature Field
│   ├── bcType.{cpp,h}               # Boundary Conditions
│   ├── fieldVal.{cpp,h}             # Helper Structure for Temperature Field
│   ├── surfaceScalarField.{cpp,h}   # Flux
│   ├── velocityField.{cpp,h}        # Fan
│   └── FileIO.{cpp,h}               # Helper Structure for Reading Inputs
├── input/                    # Configuration files
│   ├── C.inp                        # Physical Properties
│   ├── mesh.inp                     # Spatial Parameters
│   ├── nonUniformGrid.inp           # Grid Spacing Parameters
│   └── velocity.inp                 # Velocity Inputs
├── output/                   # Results and Post-Processing Scripts
│   └── SurfacePlots.m               # Animates Concentration in Membrane
├── makefile                 # Build Configuration
└──README.md                 # This File
```
