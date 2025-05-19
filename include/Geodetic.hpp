//$Header$
//------------------------------------------------------------------------------
// Geodetic
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Geodetic.hpp
 * @brief Archivo cabecera de operacion Geodetic.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _Geodetic_
#define _Geodetic_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// tuple<double, double, double> Geodetic(Matrix& r)
//------------------------------------------------------------------------------
/**
 * @brief Operacion Geodetic.
 *
 * @param [in] r Vector columna posicion en metros.
 * @return Coordenadas geodésicas (longitud [rad], latitud [rad], altitud [m]) a partir del vector de posición dado (r [m]).
 */
//------------------------------------------------------------------------------
tuple<double, double, double> Geodetic(Matrix& r);

#endif