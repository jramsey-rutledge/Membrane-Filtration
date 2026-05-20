// header file

// header guards
#ifndef NONUNIFORMGRID_H
#define NONUNIFORMGRID_H

#include <fstream>
#include <sstream>
#include <vector>
#include "FileIO.h"
#include "meshType.h"


// structure definitons
class nonUniformGrid : public meshType {
    public:
    double beta_x, beta_y;        // clustering strength

    void create() override;

};

#endif