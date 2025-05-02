//$Header$
//------------------------------------------------------------------------------
// PrecMatrix
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file PrecMatrix.hpp
 * @brief Archivo cabecera de operacion PrecMatrix.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _PrecMatrix_
#define _PrecMatrix_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& PrecMatrix(double Mjd_1, double Mjd_2)
//------------------------------------------------------------------------------
/**
 * @brief Operacion PrecMatrix.
 *
 * @param [in] Mjd_1 Epoch dada (Fecha TT juliana modificada).
 * @param [in] Mjd_2 Epoch a procesar (Fecha TT juliana modificada).
 * @return Matriz de transformacion de Precesion.
 */
//------------------------------------------------------------------------------
Matrix& PrecMatrix(double Mjd_1, double Mjd_2);

#endif