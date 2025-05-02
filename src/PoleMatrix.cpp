//$Source$
//------------------------------------------------------------------------------
// PoleMatrix
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file PoleMatrix.cpp
* @brief Programacion de operacion PoleMatrix.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\PoleMatrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"

Matrix& PoleMatrix(double xp, double yp){

    return R_y(-xp) * R_x(-yp);
}
