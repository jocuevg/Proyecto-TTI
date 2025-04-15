//$Header$
//------------------------------------------------------------------------------
// EccAnom
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file EccAnom.hpp
 * @brief Archivo cabecera de operacion EccAnom.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _EccAnom_
#define _EccAnom_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// double EccAnom (double M, double e)
//------------------------------------------------------------------------------
/**
 * @brief Operacion EccAnom.
 *
 * @param [in] M Anomalía media en radianes.  
 * @param [in] E Excentricidad de la órbita [0,1].
 * @return Anomalía excéntrica en radianes.
 */
//------------------------------------------------------------------------------
double EccAnom (double M, double e);

#endif