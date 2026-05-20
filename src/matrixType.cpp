#include "matrixType.h"
#include "meshType.h"
#include "surfaceScalarField.h"
#include "volScalarField.h"
#include "velocityField.h"
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include "Eigen/Sparse"
#include "Eigen/SparseLU"
using namespace std;

// matrixType Member Functions Defined 
// Initialize
matrixType::matrixType(const meshType & Mesh){
    Nx = Mesh.Nx ;
    Ny = Mesh.Ny ;
    Nx2 = Nx + 2;
    Ny2 = Ny + 2;
    B.assign ( Nx2 * Ny2 , 0.0);
}

// add a matrix entry at row , col
void matrixType::addEntry ( int row , int col , double val ){
    ROW.push_back ( row );
    COL.push_back ( col );
    VAL.push_back ( val );
}

// set RHS at cell (i , j )
void matrixType::setRHS ( int i , int j , double val ){
    B[ i + j * Nx2 ] = val ;
}

// get RHS at cell (i , j )
double matrixType::getRHS ( int i , int j ) const {
    return B [ i + j * Nx2 ];
}

// Function to build a sparse matrix
void matrixType::BuildA(const meshType& Mesh, const surfaceScalarField& flux, const volScalarField& T) {
    ROW.clear();
    COL.clear();
    VAL.clear();

    int row ;
    double a,b,g;

    // south / north boundary conditions
    for (int i = 1; i <= Nx; i++) {
        // south
        row = i ;
        T.BC.south(Mesh.xc[i],0.0,a,b,g);
        addEntry(row,row,0.5*a+b/Mesh.dy[0]);
        addEntry(row,row + Nx2,0.5*a-b/Mesh.dy[0]);

        // north
        row = (Ny2-1)* Nx2 + i;
        T.BC.north(Mesh.xc[i],0.0,a,b,g);
        addEntry(row,row,0.5*a+b/Mesh.dy[Ny-1]);
        addEntry(row,row - Nx2,0.5*a-b/Mesh.dy[Ny-1]);
    }

    // east / west boundary conditions
    for (int j = 1; j <= Ny; j++) {
    // west 
        row = j * Nx2;
        T.BC.west(Mesh.yc[j],0.0,a,b,g);
        addEntry(row,row,0.5*a+b/Mesh.dx[0]);
        addEntry(row,row+1, 0.5*a-b/Mesh.dx[0]);

    // east
        row = (j+1) * Nx2 - 1;
        //T.BC.east(Mesh.yc[j],0.0,a,b,g);
        //addEntry(row,row,a*0.5+b/Mesh.dx[Nx-1]);
        //addEntry(row,row-1,a*0.5-b/Mesh.dx[Nx-1]);
        addEntry(row,row, 1.0);
        addEntry(row,row-1,-2.0);
        addEntry(row,row-2,1.0);
    }

    // 5-point stencil
    for (int i=1; i<=Nx; i++) {
        for (int j=1; j<=Ny; j++) {
            row = j*Nx2 + i;
            
            // Geometry
            double dxe = Mesh.xc[i+1] - Mesh.xc[i];
            double dxw = Mesh.xc[i]   - Mesh.xc[i-1];
            double dyn = Mesh.yc[j+1] - Mesh.yc[j];
            double dys = Mesh.yc[j]   - Mesh.yc[j-1];

            double Ae = Mesh.yf[j] - Mesh.yf[j-1];
            double Aw = Ae;
            double An = Mesh.xf[i] - Mesh.xf[i-1];
            double As = An;

            // Diffusion
            double ae = T.k * Ae / dxe;
            double aw = T.k * Aw / dxw;
            double an = T.k * An / dyn;
            double as = T.k * As / dys;

            // Fluxes
            double phi_e = flux.gete(i,j);
            double phi_w = flux.getw(i,j);
            double phi_s = flux.gets(i,j);
            double phi_n = flux.getn(i,j);

            double Fe = T.cp * phi_e * Ae;
            double Fw = T.cp * phi_w * Aw;
            double Fn = T.cp * phi_n * An;
            double Fs = T.cp * phi_s * As;

            // Central differencing convection
            double cE =  0.5*Fe;
            double cW = -0.5*Fw;
            double cN =  0.5*Fn;
            double cS = -0.5*Fs;
            double cP = 0.5*(Fe-Fw+Fn-Fs);

            double ap = ae+aw+an+as;

            // Matrix
            addEntry(row,row,ap+cP);
            addEntry(row,row+1,-ae+cE);
            addEntry(row,row-1,-aw+cW);
            addEntry(row,row+Nx2,-an+cN);
            addEntry(row,row-Nx2,-as+cS);
        }
    }

    // extraneous nodes
    addEntry(0,0,1.0);
    addEntry(Nx+1,Nx+1,1.0);
    addEntry((Ny2-1)*Nx2,(Ny2-1)*Nx2,1.0);
    addEntry(Nx2*Ny2-1,Nx2*Ny2-1,1.0);

}

// Build right - hand side vector
void matrixType::BuildRHS( const meshType& Mesh, const volScalarField& T, const velocityField& U ){
    std::fill(B.begin(), B.end(), 0.0);
    double a,b,g;
    
    // south / north boundary conditions
    for ( int i = 1; i <= Nx ; i++){
    
        T.BC.south(Mesh.xc[i],0.0,a,b,g);
        setRHS(i,0,g);

        T.BC.north(Mesh.xc[i],0.0,a,b,g);
        setRHS(i,Ny2-1,g);
        
    }

    // east / west boundary conditions
    for ( int j = 1; j <= Ny ; j ++){

        T.BC.east(Mesh.yc[j],0.0,a,b,g);
        setRHS(Nx2-1,j,g);

        T.BC.west(Mesh.yc[j],0.0,a,b,g);
        setRHS(0,j,g);
    }

    // internal nodes
    for (int i = 1; i <= Nx; i++){
        for (int j = 1; j <= Ny; j++){
            setRHS(i, j, 0.0);
        }
    }

    // unusued parameters
    double Vm = U.Vm;
    Vm = Vm;
}

// Solve Matrix
void matrixType::solve(volScalarField& T){
  int Ntot = Nx2 * Ny2;

  // load matrix
  Eigen::SparseMatrix<double> A(Ntot, Ntot);
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(VAL.size());
  for (size_t k = 0; k < VAL.size(); ++k){
    triplets.emplace_back(ROW[k], COL[k], VAL[k]);
  }
  A.setFromTriplets(triplets.begin(), triplets.end());

  // load right-hand-side
  Eigen::VectorXd RHS(Ntot);
  for (int i = 0; i < Ntot; ++i)
  {
    RHS[i] = B[i];
  }

  // compute the LU factorization of A
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  if(solver.info() != Eigen::Success) 
  {
    std::cerr << "Decomposition failed!\n";
    return;
  }

  // solve
  Eigen::VectorXd C = solver.solve(RHS);
  if(solver.info() != Eigen::Success) 
  {
    std::cerr << "Solving failed!\n";
    return;
  }

  // load solution into volScalarField T
  for (int j = 0; j < Ny2; ++j)
  {
    for (int i = 0; i < Nx2; ++i)
    {
      T.set(i, j, C[i + j*Nx2]);
    }
  }
}
