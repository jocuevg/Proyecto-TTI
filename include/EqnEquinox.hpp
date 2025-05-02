//$Header$
//------------------------------------------------------------------------------
// EqnEquinox
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file EqnEquinox.hpp
 * @brief Archivo cabecera de operacion EqnEquinox.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _EqnEquinox_
#define _EqnEquinox_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// double EqnEquinox(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * @brief Operacion EqnEquinox.
 *
 * @param [in] Mjd_TT Fecha juliana modificada (hora terrestre).
 * @return Ecuacion de equinocio.
 */
//------------------------------------------------------------------------------
double EqnEquinox(double Mjd_TT);

#endif