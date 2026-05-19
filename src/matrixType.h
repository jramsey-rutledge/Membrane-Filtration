// header file

// header guards
#ifndef MATRIXTYPE_H
#define MATRIXTYPE_H

// libraries
#include <vector>
class meshType;
class surfaceScalarField;
class volScalarField;
class velocityField;


// structure definitions
struct matrixType {
    int Nx , Ny ;               // number of cells in x and y
    int Nx2 , Ny2 ;             // include ghost cells  
    std::vector<int> ROW ;      // row indices of nonzeros
    std::vector<int> COL ;      // column indices of nonzeros
    std::vector<double> VAL ;   // values of nonzeros
    std::vector<double> B ;     // right - hand side vector

    // constructor
    matrixType( const meshType & Mesh );

    // add a matrix entry at row , col
    void addEntry ( int row , int col , double val );

    // set RHS at cell (i , j )
    void setRHS ( int i , int j , double val );

    // get RHS at cell (i , j )
    double getRHS ( int i , int j ) const ;

    // build sparse matrix
    void BuildA(const meshType& Mesh, const surfaceScalarField& flux, const volScalarField& T);
    
    // build RHS vector
    void BuildRHS (const meshType& Mesh, const volScalarField& T, const velocityField& U);

    // solve
    void solve(volScalarField& T);

};


#endif