//$Header$
//------------------------------------------------------------------------------
// gast
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file gast.hpp
 * @brief Archivo cabecera de operacion gast.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _gast_
#define _gast_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// double gast(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * @brief Operacion gast.
 *
 * @param [in] Mjd_UT1 Fecha UT1 juliana modificada.
 * @return GAST en radianes.
 */
//------------------------------------------------------------------------------
double gast(double Mjd_UT1);

#endif