//$Source$
//------------------------------------------------------------------------------
// GHAMatrix
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file GHAMatrix.cpp
* @brief Programacion de operacion GHAMatrix.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\GHAMatrix.hpp"
#include "..\include\R_z.hpp"
#include "..\include\gast.hpp"

Matrix& GHAMatrix(double Mjd_UT1){

    return R_z( gast(Mjd_UT1) );
}