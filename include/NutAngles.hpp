//$Header$
//------------------------------------------------------------------------------
// NutAngles
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file NutAngles.hpp
 * @brief Archivo cabecera de operacion NutAngles.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _NutAngles_
#define _NutAngles_

#include <cmath>
#include <tuple>
#include "matrix.hpp"
using namespace std;

//------------------------------------------------------------------------------
// tuple<double,double> NutAngles (double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * @brief Operacion NutAngles.
 *
 * @param [in] Mjd_TT Fecha juliana modificada (TT).  
 * @return Angulos de nutacion.
 */
//------------------------------------------------------------------------------
tuple<double,double> NutAngles (double Mjd_TT);

#endif