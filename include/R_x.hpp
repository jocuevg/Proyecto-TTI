//$Header$
//------------------------------------------------------------------------------
// R_x
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file R_x.hpp
 * @brief Archivo cabecera de operacion R_x.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _R_x_
#define _R_x_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& R_x(double angle)
//------------------------------------------------------------------------------
/**
 * @brief Operacion R_x.
 *
 * @param [in] angle angulo de rotacion.
 * @return Vector resultado.
 */
//------------------------------------------------------------------------------
Matrix& R_x(double angle);

#endif