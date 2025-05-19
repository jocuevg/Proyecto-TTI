//$Source$
//------------------------------------------------------------------------------
// Geodetic
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Geodetic.cpp
 * @brief Programacion de operacion Geodetic.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#include "..\include\Geodetic.hpp"
#include "..\include\Sat_const.hpp"

tuple<double, double, double> Geodetic(Matrix& r)
{
    double R_equ = R_Earth;
    double f     = f_Earth;

    double epsRequ = __DBL_EPSILON__*R_equ;        // Convergence criterion
    double e2      = f*(2.0-f);        // Square of eccentricity

    double X = r(1,1);                   // Cartesian coordinates
    double Y = r(2,1);
    double Z = r(3,1);
    double rho2 = X*X + Y*Y;           // Square of distance from z-axis

    double lon,lat,h;

    // Check validity of input data
    if (norm(transpose(r))==0.0){
        cout << "invalid input in Geodetic constructor\n";
		exit(EXIT_FAILURE);
        lon = 0.0;
        lat = 0.0;
        h   = -R_Earth;
    }

    // Iteration 
    double dZ = e2*Z;
    double ZdZ,Nh,SinPhi,N,dZ_new;

    while(1){
        ZdZ    =  Z + dZ;
        Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
        SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
        N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
        dZ_new =  N*e2*SinPhi;
        if ( abs(dZ-dZ_new) < epsRequ )
            break;
        
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    lon = atan2 ( Y, X );
    lat = atan2 ( ZdZ, sqrt(rho2) );
    h   = Nh - N;

    return tie(lon,lat,h);
}