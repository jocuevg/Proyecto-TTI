//$Header$
//------------------------------------------------------------------------------
// JPL_Eph_DE430
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file JPL_Eph_DE430.hpp
 * @brief Archivo cabecera de operacion JPL_Eph_DE430.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _JPL_Eph_DE430_
#define _JPL_Eph_DE430_

#include <cmath>
#include <tuple>
#include "matrix.hpp"
#include "global.hpp"
using namespace std;

//------------------------------------------------------------------------------
// tuple<Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&> JPL_Eph_DE430 (double Mjd_TDB);
//------------------------------------------------------------------------------
/**
 * @brief Operacion JPL_Eph_DE430.
 *
 * @param [in] Mjd_TDB Fecha juliana modificada de TDB. 
 * @return r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
 *  r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
 *  r_Sun(geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)).
 */
//------------------------------------------------------------------------------
tuple<Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&> JPL_Eph_DE430 (double Mjd_TDB);

#endif