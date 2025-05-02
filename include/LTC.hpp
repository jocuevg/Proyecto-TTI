//$Header$
//------------------------------------------------------------------------------
// LTC
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file LTC.hpp
 * @brief Archivo cabecera de operacion LTC.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _LTC_
#define _LTC_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& LTC(double lon, double lat)
//------------------------------------------------------------------------------
/**
 * @brief Operacion LTC.
 *
 * @param [in] lon Longitud.
 * @param [in] lat Latitud.
 * @return Matriz de rotacion terrestre.
 */
//------------------------------------------------------------------------------
Matrix& LTC(double lon, double lat);

#endif