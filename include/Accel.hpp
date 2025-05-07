//$Header$
//------------------------------------------------------------------------------
// Accel
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Accel.hpp
 * @brief Archivo cabecera de operacion Accel.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _Accel_
#define _Accel_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& Accel(double x,Matrix& Y)
//------------------------------------------------------------------------------
/**
 * @brief Operacion Accel.
 *
 * @param [in] x Hora terrestre.
 * @param [in] Y Vector columna del estado del satelite en el sistema ICRF/EME2000.
 * @return Vector (fila por conveniencia) aceleracion (a=d^2r/dt^2) en el sistema ICRF/EME2000.
 */
//------------------------------------------------------------------------------
Matrix& Accel(double x,Matrix& Y);

#endif