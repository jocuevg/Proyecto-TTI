//$Header$
//------------------------------------------------------------------------------
// AccelPointMass
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file AccelPointMass.hpp
 * @brief Archivo cabecera de operacion AccelPointMass.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _AccelPointMass_
#define _AccelPointMass_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// double AccelPointMass(Matrix r, Matrix s,double GM)
//------------------------------------------------------------------------------
/**
 * @brief Operacion AccelPointMass.
 *
 * @param [in] r Satellite position vector.
 * @param [in] s Point mass position vector.
 * @param [in] GM Gravitational coefficient of point mass.
 * @return aceleracion.
 */
//------------------------------------------------------------------------------
Matrix AccelPointMass(Matrix r, Matrix s,double GM);

#endif