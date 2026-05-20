// header file

// header guards
#ifndef UNIFORMGRID_H
#define UNIFORMGRID_H

#include <fstream>
#include <sstream>
#include <vector>
#include "FileIO.h"
#include "meshType.h"


// structure definitons
class uniformGrid : public meshType {
    public:

// function prototypes
/**
 * @brief constructor of Mesh grid
 * 
 * @param Mesh passed by reference to make a copy
 */
void create() override;




};

#endif