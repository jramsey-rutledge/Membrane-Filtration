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

    // Assign x values
    double dx_uniform = Lx / Nx;
    dx.assign(Nx, dx_uniform);

    xf.assign(Nx+1,0.0);
    for(int i=0; i<=Nx; i++){ xf[i] = i * dx_uniform; }

    xc.assign(Nx+2,0.0);
    for(int i=1; i<=Nx; i++){ xc[i] = 0.5 * (xf[i-1] + xf[i]); }

    // ghost cells
    xc[0] = xc[1] - dx_uniform;
    xc[Nx+1] = xc[Nx] + dx_uniform;

    // Assign y values
    double dy_uniform = Ly / Ny;
    dy.assign(Ny, dy_uniform);

    yf.assign(Ny+1,0.0);
    for(int j=0; j<=Ny; j++){ yf[j] = j * dy_uniform; }

    yc.assign(Ny+2,0.0);
    for(int j=1; j<=Ny; j++){ yc[j] = 0.5 * (yf[j-1] + yf[j]); }

    // ghost cells
    yc[0] = yc[1] - dy_uniform;
    yc[Ny+1] = yc[Ny] + dy_uniform;

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

