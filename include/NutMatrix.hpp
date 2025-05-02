//$Header$
//------------------------------------------------------------------------------
// NutMatrix
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file NutMatrix.hpp
 * @brief Archivo cabecera de operacion NutMatrix.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _NutMatrix_
#define _NutMatrix_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& NutMatrix(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * @brief Operacion NutMatrix.
 *
 * @param [in] Mjd_TT Fecha juliana modificada (hora terrestre).
 * @return Matriz de nutacion.
 */
//------------------------------------------------------------------------------
Matrix& NutMatrix(double Mjd_TT);

#endif