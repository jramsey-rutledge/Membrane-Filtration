// header guards
#ifndef VELOCITYFIELD_H
#define VELOCITYFIELD_H

#include "FileIO.h"
#include "meshType.h"

struct velocityField{
// --- member data --- //
double Uin ; // mean inlet velocity
double Vm ; // transmembrane velocity
double h; // channel half - height

// --- member functions --- //
velocityField ( const meshType & Mesh );
double compU ( double x , double y ) const ;
double compV ( double x , double y ) const ;

};
#endif