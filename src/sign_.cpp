//$Source$
//------------------------------------------------------------------------------
// sign_
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file sign_.cpp
 * @brief Programacion de operacion sign_.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#include "..\include\sign_.hpp"

double sign_(double a, double b)
{
    if (b >= 0.0)
        return fabs(a);
    else
        return -fabs(a);
}
