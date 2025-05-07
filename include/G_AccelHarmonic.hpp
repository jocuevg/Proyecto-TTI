//$Header$
//------------------------------------------------------------------------------
// G_AccelHarmonic
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file G_AccelHarmonic.hpp
 * @brief Archivo cabecera de operacion G_AccelHarmonic.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _G_AccelHarmonic_
#define _G_AccelHarmonic_

#include <cmath>
#include "matrix.hpp"
#include "global.hpp"

//------------------------------------------------------------------------------
// Matrix& G_AccelHarmonic(Matrix& r,Matrix& U,int n_max,int m_max )
//------------------------------------------------------------------------------
/**
 * @brief Operacion G_AccelHarmonic.
 *
 * @param [in] r Vector columna posicion del satelite en el sistema true-of-date.
 * @param [in] U Matriz de transformación al sistema body-fixed.
 * @param [in] n_max Grado del modelo de gravedad.
 * @param [in] m_max Orden del modelo de gravedad.
 * @return Gradiente (G=da/dr) en el sistema true-of-date.
 */
//------------------------------------------------------------------------------
Matrix& G_AccelHarmonic(Matrix& r,Matrix& U,int n_max,int m_max );

#endif