/* MEGN 651: Modeling Membrane Filtration
 *
 * Author: Jenna Ramsey-Rutledge
 *
 * In this code, we build a computational model of concentration polarization and use it to explore how the phenomenon varies across filtration regimes. 
 */

// Libraries
#include "meshType.h"
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
    meshType Mesh;

    // Create Matrix A and RHS vector
    volScalarField C("C", Mesh);

    // Create velocity field
    velocityField U(Mesh);
    C.BC.Vm = U.Vm;

    // Create flux
    surfaceScalarField flux;
    flux.allocate(Mesh);
    flux.computePhix(Mesh,C,U);
    flux.computePhiy(Mesh,C,U);

    matrixType Matrix(Mesh);          // initialize matrix
    Matrix.BuildA(Mesh,flux,C);       // build matrix
    Matrix.BuildRHS(Mesh,C,U);        // build right hand side
    Matrix.solve(C);                  // solve using Eigen

    C.Save(); // save field to file

    // // Debugging for manufactured solution
    // double error = 0.0;

    // for (int i = 1; i <= Mesh.Nx; i++){
    //     for (int j = 1; j <= Mesh.Ny; j++){
    //         double x = Mesh.xc[i];
    //         double y = Mesh.yc[j];

    //         double exact = sin(x)*sin(y);
    //         double diff  = C.get(i,j) - exact;

    //         error += diff*diff;
    //     }
    // }

    // error = sqrt(error / (Mesh.Nx*Mesh.Ny));
    // cout << "L2 Error = " << error << endl;
return 0;
}