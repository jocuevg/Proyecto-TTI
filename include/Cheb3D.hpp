//$Header$
//------------------------------------------------------------------------------
// Cheb3D
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Cheb3D.hpp
 * @brief Archivo cabecera de operacion Cheb3D.
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
 Matrix Cheb3D(double t,double N,double Ta,double Tb,Matrix Cx,Matrix Cy,Matrix Cz);

#endif