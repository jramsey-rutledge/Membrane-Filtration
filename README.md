# Membrane-Filtration Model
This code solves the steady advection-diffusion equation for a membrane filtration system. 

**Author:**
- Jenna Ramsey-Rutledge

## Background

## Physics and Mathematical Model

### Governing Equations
The concentration $C(x,y)$ satisfies the steady advection-diffusion equation written in integral form, 

$$\int_{CS} (\rho c_{p} \vec{u} C) \cdot \vec{n} dA = \int_{CS} k (\nabla C \cdot \vec{n}) dA + \int_{CV} f dV,$$

where $\rho$, $c_p$, and $k$ are the density, specific heat, and theraml conductivity of air, which we assume are constants. The velocity field is defined $\vec{u}=(u,v)$. The source term $f$ is included both for verification using the Method of Manufactured Solutions (MMS) and for modeling internal heat generation such as a heater.


### Initial and Boundary Conditions
This solver supports both Dirichlet and Neumann boundary conditions through a generalized Robin boundary-condition framework,

$$aT + b\frac{\partial T}{\partial n} = g.$$


## Numerical Methods 
The governing equations are discretized using the finite-volume method on both a uniform and non-uniform grid with cell-centered unknowns. The spatial discretization utilizes central differencing. The diffusion terms are handled with second-order central differences and the advection terms use a centered finite-volume flux formulation. The face fluxes are computed via linear interpolation of the velocity values on cell faces. 



## Verification



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
Unsteady-Fan-MOdel/
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
