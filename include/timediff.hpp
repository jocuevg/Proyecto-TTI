//$Header$
//------------------------------------------------------------------------------
// timediff
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file timediff.hpp
 * @brief Archivo cabecera de operacion timediff.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _timediff_
#define _timediff_

#include <cmath>
#include <tuple>
using namespace std;

//------------------------------------------------------------------------------
// tuple<double,double,double,double,double> timediff (double UT1_UTC, double TAI_UTC)
//------------------------------------------------------------------------------
/**
 * @brief Operacion timediff.
 *
 * @param [in] UT1_UTC Primer tiempo. 
 * @param [in] TAI_UTC Segundo tiempo.  
 * @return Diferencias de tiempo.
 */
//------------------------------------------------------------------------------
tuple<double,double,double,double,double> timediff (double UT1_UTC, double TAI_UTC);

#endif