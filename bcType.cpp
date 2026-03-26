#include "bcType.h"
#include <fstream>
#include <iostream>

// --- constructor --- //
bcType :: bcType (){
    // read Vm from velocity . inp
    std::ifstream file1;
    file1.open ("velocity.inp");
    Vm = readDouble ( file1 );
    file1.close ();

    // read k from file C. inp
    std::ifstream file2;
    file2.open("C.inp");
    readDouble(file2); readDouble(file2); // disregard first two entries
    k = readDouble(file2);
    file2.close();

    // display for verification
    std :: cout << "- -- BC data -- -" << std :: endl ;
    std :: cout << " Vm = " << Vm << std :: endl ;
    std :: cout << "k = " << k << std :: endl ;

}

// west boundary condition
void bcType::west( double y, double t, double& a, double& b, double& g) const {
    a = 1.0;
    b = 0.0;
    g = 1.0;
    t=t; y=y;

}

// east boundary condition
void bcType::east( double y, double t, double& a, double& b, double& g) const {
    a = 0.0; // coefficient multiplying T
    b = 1.0; // coefficient multiplying dT / dx
    g = 0.0; // prescribed value
    t=t;
    y=y;
}

// north boundary condition
void bcType::north( double x, double t, double& a, double& b, double& g) const {
    a = -Vm;
    b = k;
    g = 0.0;
    x=x;
    t=t;
}

// south boundary condition
void bcType::south( double x, double t, double& a, double& b, double& g) const {
    a = Vm;
    b = -k;
    g = 0.0;
    x=x;
    t=t;
}