//$Header$
//------------------------------------------------------------------------------
// MeasUpdate
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file MeasUpdate.hpp
 * @brief Archivo cabecera de operacion MeasUpdate.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _MeasUpdate_
#define _MeasUpdate_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// tuple<Matrix&,Matrix&,Matrix&> MeasUpdate(Matrix& x,double z,double g,double s,Matrix& G,Matrix& P,int n)
//------------------------------------------------------------------------------
/**
 * @brief Operacion MeasUpdate.
 *
 * @param [in] x Vector columna.
 * @param [in] z Double.
 * @param [in] g Double.
 * @param [in] s Double.
 * @param [in] G Vector fila.
 * @param [in] P Matrix.
 * @param [in] n Entero.
 * @return Tupla de matrices con K (vector columna), x y P.
 */
//------------------------------------------------------------------------------
tuple<Matrix&,Matrix&,Matrix&> MeasUpdate(Matrix& x,double z,double g,double s,Matrix& G,Matrix& P,int n);

#endif