/* MEGN 651: Modeling Membrane Filtration
 *
 * Author: Jenna Ramsey-Rutledge
 *
 * In this code, we build a computational model of concentration polarization and use it to explore how the phenomenon varies across filtration regimes. 
 */

// Libraries
#include "nonUniformGrid.h"
#include "uniformGrid.h"
#include "matrixType.h"
#include "volScalarField.h"
#include "surfaceScalarField.h"
#include "velocityField.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath> 
#include <vector>

using namespace std;


// Main Function
int main() {

    // Create Mesh
    //uniformGrid Mesh;      // uniform mesh   
    nonUniformGrid Mesh;   // non uniform mesh
    Mesh.create(); 

    // Create Matrix A and RHS vector
    volScalarField C("C", Mesh);

    // Create velocity field
    velocityField U(Mesh);
    C.BC.Vm = U.Vm;

    // Create flux
    surfaceScalarField flux;
    flux.allocate(Mesh);
    flux.computePhix(Mesh,U);
    flux.computePhiy(Mesh,U);

    matrixType Matrix(Mesh);          // initialize matrix
    Matrix.BuildA(Mesh,flux,C);       // build matrix
    Matrix.BuildRHS(Mesh,C,U);        // build right hand side
    Matrix.solve(C);                  // solve using Eigen

    C.Save(); // save field to file

return 0;
}