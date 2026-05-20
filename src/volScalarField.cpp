// Libraries
#include "meshType.h"
#include "volScalarField.h"
#include "FileIO.h"
#include "bcType.h"
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <math.h>

using namespace std;

void volScalarField :: read(){
    std :: ifstream file ;
    file.open ("input/" + name + ".inp");
    rho = readDouble ( file );
    cp = readDouble ( file );
    k = readDouble ( file );
    file.close ();
}

void volScalarField::dispInputs() const {
    std::cout << "- -- " << name << " inputs - --" << "\n";
    std::cout << " rho = " << rho << "\n ";
    std::cout << " cp = " << cp << "\n ";
    std::cout << "k = " << k << "\n ";
}

void volScalarField::allocate(const meshType& Mesh) {
    Nx = Mesh.Nx;
    Ny = Mesh.Ny;
    Ntot = (Nx+2)*(Ny+2);
    values.assign(Ntot, 0.0);
}

// Constructors
volScalarField::volScalarField(){
    name = "T";
    read();
    dispInputs();
}

volScalarField::volScalarField(string fieldName){
    name = fieldName;
    read();
    dispInputs();
}

volScalarField::volScalarField(string fieldName, const meshType& Mesh){
    name = fieldName;
    read();
    dispInputs();
    allocate(Mesh);
}

// Getters and Setters
double volScalarField::get(int i, int j) const {
    return values[i+j*(Nx+2)];
}

void volScalarField::set(int i, int j, double val) {
    values[i+j*(Nx+2)] = val;
}

// Save final results
void volScalarField::Save() {
    std::string filename = "output/" + name + "_" + std::to_string(Nx) + ".dat";
    std::ofstream cfile(filename);
    cfile << std::scientific << std::setprecision(16);
    for (int i = 0; i < Ntot; i++) {cfile << values[i] << "\n";}
    cfile.close();
}