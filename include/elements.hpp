//$Header$
//------------------------------------------------------------------------------
// elements
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file elements.hpp
 * @brief Archivo cabecera de operacion elements.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _elements_
#define _elements_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// tuple<double,double,double,double,double,double,double> elements(Matrix& y)
//------------------------------------------------------------------------------
/**
 * @brief Operacion elements.
 *
 * @param [in] vec1 Vector de estado (x,y,z,vx,vy,vz).
 * @return p - semilatus rectum [m], a - Semimajor axis, e - Eccentricity, i - Inclination [rad], 
 * Omega - Longitude of the ascending node [rad], omega - Argument of pericenter [rad], M - mean anomaly [rad].
 */
//------------------------------------------------------------------------------
tuple<double,double,double,double,double,double,double> elements(Matrix& y);

#endif