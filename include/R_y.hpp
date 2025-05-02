//$Header$
//------------------------------------------------------------------------------
// R_y
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file R_y.hpp
 * @brief Archivo cabecera de operacion R_y.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _R_y_
#define _R_y_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& R_y(double angle)
//------------------------------------------------------------------------------
/**
 * @brief Operacion R_y.
 *
 * @param [in] angle angulo de rotacion.
 * @return Vector resultado.
 */
//------------------------------------------------------------------------------
Matrix& R_y(double angle);

#endif