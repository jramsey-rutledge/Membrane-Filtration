#include "uniformGrid.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <math.h>
using namespace std;

void uniformGrid::create(){
    // Read inputs
    ifstream infile;
    infile.open("input/mesh.inp");

    if (!infile) {
        cout << "Error opening File\n";
        exit(1);
    }

    Nx = readInt(infile);
    Ny = readInt(infile);
    Lx = readDouble(infile);
    Ly = readDouble(infile); 

    // Grid Spacing
    dx.assign(Nx,Lx/Nx);
    dy.assign(Ny,Ly/Ny);


    // Display inputs
    cout << "- -- Mesh inputs - --" << endl ;
    cout << "Lx = " << Lx << endl;
    cout << "Ly = " << Ly << endl;
    cout << "Nx = " << Nx << endl;
    cout << "Ny = " << Ny << endl;

    // Assign values
    xc.assign(Nx+2,0.0);
    for (int i=0; i<Nx+2;i++){
        xc[i] = -0.5*dx[i] + i*dx[i];
    }

    yc.assign(Ny+2,0.0);
    for (int i=0; i<Ny+2;i++){
        yc[i] = -0.5*dy[i] + i*dy[i];
    }

    xf.assign(Nx+1,0.0);
    for (int k=0; k<Nx+1; ++k){
        xf[k] = k * dx[k];
    }

    yf.assign(Ny+1,0.0);
    for (int k=0; k<Ny+1; ++k){
        yf[k] = k * dy[k];
    }

    // write values
    std::ofstream xfile("output/x.dat");
    for (int i = 0; i < Nx+2; i++){
        xfile << xc[i] << "\n";
    }
    xfile.close();

    std::ofstream yfile("output/y.dat");
    for (int i = 0; i < Ny+2; i++) {
        yfile << yc[i] << "\n";
    }
    yfile.close();
}

