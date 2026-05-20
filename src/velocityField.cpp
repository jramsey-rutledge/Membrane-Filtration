# include "velocityField.h"
#include <iostream>

velocityField :: velocityField ( const meshType & Mesh ){
    h = 0.5* Mesh.Ly ;
    std :: ifstream file ;
    file . open ("input/velocity.inp");
    Uin = readDouble ( file );
    Vm = readDouble ( file );
    file.close ();

    std :: cout << "- -- velocity inputs - --" << std :: endl ;
    std :: cout << " Uin = " << Uin << std :: endl ;
    std :: cout << " Vm = " << Vm << std :: endl ;
}

double velocityField :: compU ( double x , double y) const {
    double xhat = x/h;
    double yhat = ( y - h )/ h;
    return 1.5* Uin *(1.0 - Vm * xhat / Uin )*(1.0 - yhat * yhat );
}

double velocityField :: compV ( double x , double y) const {
    double yhat = ( y - h )/ h;
    return 0.5* Vm * yhat *(3.0 - yhat * yhat );
    x=x;
}