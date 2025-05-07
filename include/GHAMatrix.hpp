//$Header$
//------------------------------------------------------------------------------
// GHAMatrix
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file GHAMatrix.hpp
 * @brief Archivo cabecera de operacion GHAMatrix.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _GHAMatrix_
#define _GHAMatrix_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& GHAMatrix(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * @brief Operacion GHAMatrix.
 *
 * @param [in] Mjd_UT1 Fecha UT1 juliana modificada.
 * @return Matriz angulo de la hora de Greenwich.
 */
//------------------------------------------------------------------------------
Matrix& GHAMatrix(double Mjd_UT1);

#endif