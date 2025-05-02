//$Source$
//------------------------------------------------------------------------------
// EqnEquinox
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file EqnEquinox.cpp
* @brief Programacion de operacion EqnEquinox.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\EqnEquinox.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\MeanObliquity.hpp"

double EqnEquinox(double Mjd_TT){
    // Nutation in longitude and obliquity
    auto[dpsi, deps] = NutAngles (Mjd_TT);

    // Equation of the equinoxes
    return dpsi * cos ( MeanObliquity(Mjd_TT) );
}
