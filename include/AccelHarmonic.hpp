//$Header$
//------------------------------------------------------------------------------
// AccelHarmonic
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file AccelHarmonic.hpp
 * @brief Archivo cabecera de operacion AccelHarmonic.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _AccelHarmonic_
#define _AccelHarmonic_

#include <cmath>
#include "matrix.hpp"
#include "global.hpp"

//------------------------------------------------------------------------------
// Matrix& AccelHarmonic(Matrix& r, Matrix& E,int n_max,int m_max)
//------------------------------------------------------------------------------
/**
 * @brief Operacion AccelHarmonic.
 *
 * @param [in] r Vector posicion del satelite.
 * @param [in] E Matriz de transformación.
 * @param [in] n_max Grado maximo.
 * @param [in] m_max Orden maximo.
 * @return Aceleracion.
 */
//------------------------------------------------------------------------------
Matrix& AccelHarmonic(Matrix& r, Matrix& E,int n_max,int m_max);

#endif