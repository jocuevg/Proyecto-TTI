//$Header$
//------------------------------------------------------------------------------
// MeanObliquity
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file MeanObliquity.hpp
 * @brief Archivo cabecera de operacion MeanObliquity.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _MeanObliquity_
#define _MeanObliquity_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// double MeanObliquity(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * @brief Operacion MeanObliquity.
 *
 * @param [in] Mjd_TT Fecha juliana modificada(Tiempo terrestre).  
 * @return Oblicuidad media de la eclíptica [rad].
 */
//------------------------------------------------------------------------------
double MeanObliquity(double Mjd_TT);

#endif