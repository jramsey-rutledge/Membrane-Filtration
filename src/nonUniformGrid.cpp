#include "nonUniformGrid.h"
#include "FileIO.h"
#include "meshType.h"

#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

void nonUniformGrid::create() {
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

    // Assign vectors
    xc.assign(Nx+2,0.0);
    yc.assign(Ny+2,0.0);

    xf.assign(Nx+1,0.0);
    yf.assign(Ny+1,0.0);

    dx.assign(Nx,0.0);
    dy.assign(Ny,0.002/64); //uniform grid spacing in y-direction

    // Read Important Parameters
    std::ifstream file("input/nonUniformGrid.inp");

    if (!file) {
        std::cout << "Error Opening nonUniformGrid File\n";
        exit(1);
    }

    beta_x = readDouble(file);
    beta_y = readDouble(file);

    // Assign y faces 
    for(int i=0; i<=Ny; i++){
        double s = (double)i / Ny;
        yf[i] = 0.5 * Ly * (1.0 + tanh(beta_y*(2.0*s - 1.0)) / tanh(beta_y));
    }

    // Define cell length
    for(int i=0; i<Ny; i++){ dy[i] = yf[i+1] - yf[i]; }

    // Assign y centers
    for(int i=1; i<=Ny; i++){ yc[i] = 0.5 * (yf[i-1] + yf[i]); }

    // Deal with ghost cells
    yc[0] = yc[1] - dy[0];
   yc[Ny+1] = yc[Ny] + (yc[Ny]-yc[Ny-1]);

    // Assign uniform x faces
    for(int j=0; j<=Nx; j++){ xf[j] = j * dx[j]; }

    // Assign uniform x centers
    for(int j=0; j<=Nx+1; j++){ xc[j] = -0.5*dx[j] + j*dx[j]; }
    
    // Write x values to data file
    std::ofstream xfile("output/x.dat");
    for (int i = 0; i < Nx+2; i++){
        xfile << xc[i] << "\n";
    }
    xfile.close();

    // Write y values to data file
    std::ofstream yfile("output/y.dat");
    for (int i = 0; i < Ny+2; i++) {
        yfile << yc[i] << "\n";
    }
    yfile.close();

}