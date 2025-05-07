//$Header$
//------------------------------------------------------------------------------
// R_z
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file R_z.hpp
 * @brief Archivo cabecera de operacion R_z.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _R_z_
#define _R_z_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& R_z(double angle)
//------------------------------------------------------------------------------
/**
 * @brief Operacion R_z.
 *
 * @param [in] angle angulo de rotacion.
 * @return Matriz resultado.
 */
//------------------------------------------------------------------------------
Matrix& R_z(double angle);

#endif