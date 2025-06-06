//$Source$
//------------------------------------------------------------------------------
// gast
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file gast.cpp
* @brief Programacion de operacion gast.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\gast.hpp"
#include "..\include\gmst.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\Sat_const.hpp"

double gast(double Mjd_UT1){

    double gstime= fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), pi2 );
    if(gstime<0) gstime+=pi2;
    return gstime;
}
