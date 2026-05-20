// header file

// header guards
#ifndef MESHTYPE_H
#define MESHTYPE_H

#include <fstream>
#include <sstream>
#include <vector>
#include "FileIO.h"

// parent structure definitons
class meshType {
public:
    int Nx, Ny;                   // number of cells
    double Lx, Ly;                // size of domain
    std::vector<double> dx, dy;   // cell size
    std::vector<double> xc, yc;   // cell centers
    std::vector<double> xf, yf;   // cell faces

    virtual void create() = 0;

    virtual ~meshType() {}
};

#endif