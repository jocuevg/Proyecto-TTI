//$Header$
//------------------------------------------------------------------------------
// angl
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file angl.hpp
 * @brief Archivo cabecera de operacion angl.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _angl_
#define _angl_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// double angl(Matrix& vec1,Matrix& vec2)
//------------------------------------------------------------------------------
/**
 * @brief Operacion angl.
 *
 * @param [in] vec1 Vector columna 1.
 * @param [in] vec2 Vector columna 2.
 * @return Angulo entre los dos vectores.
 */
//------------------------------------------------------------------------------
double angl(Matrix& vec1,Matrix& vec2);

#endif