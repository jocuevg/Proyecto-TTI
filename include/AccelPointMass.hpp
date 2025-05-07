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
// Matrix& AccelPointMass(Matrix& r, Matrix& s,double GM)
//------------------------------------------------------------------------------
/**
 * @brief Operacion AccelPointMass.
 *
 * @param [in] r Vector columna posicion del satelite.
 * @param [in] s Vector columna posicion del punto de masa.
 * @param [in] GM Coeficiente gravitacional de la masa.
 * @return Vector columna aceleracion.
 */
//------------------------------------------------------------------------------
Matrix& AccelPointMass(Matrix& r, Matrix& s,double GM);

#endif