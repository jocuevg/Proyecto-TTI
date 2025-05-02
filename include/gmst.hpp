//$Header$
//------------------------------------------------------------------------------
// gmst
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file gmst.hpp
 * @brief Archivo cabecera de operacion gmst.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _gmst_
#define _gmst_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// double gmst(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * @brief Operacion gmst.
 *
 * @param [in] Mjd_UT1 Fecha UT1 juliana modificada.
 * @return GMST en radianes.
 */
//------------------------------------------------------------------------------
double gmst(double Mjd_UT1);

#endif