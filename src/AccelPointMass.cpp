//$Source$
//------------------------------------------------------------------------------
// AccelPointMass
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file AccelPointMass.cpp
* @brief Programacion de operacion AccelPointMass.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\AccelPointMass.hpp"

Matrix& AccelPointMass(Matrix& r, Matrix& s,double GM){
    Matrix& d = r - s;

    return ( d/(pow(norm(transpose(d)),3)) + s/(pow(norm(transpose(s)),3)) ) * -GM;
}
