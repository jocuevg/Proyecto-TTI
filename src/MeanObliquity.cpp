//$Source$
//------------------------------------------------------------------------------
// MeanObliquity
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file MeanObliquity.cpp
 * @brief Programacion de operacion MeanObliquity.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#include "..\include\MeanObliquity.hpp"
#include "..\include\Sat_const.hpp"

double MeanObliquity(double Mjd_TT)
{
    double T = (Mjd_TT-MJD_J2000)/36525;

    return Rad *( 84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600 );
}
