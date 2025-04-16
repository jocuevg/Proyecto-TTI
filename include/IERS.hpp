//$Header$
//------------------------------------------------------------------------------
// IERS
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file IERS.hpp
 * @brief Archivo cabecera de operacion IERS.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _IERS_
#define _IERS_

#include <cmath>
#include <tuple>
#include "matrix.hpp"
using namespace std;

//------------------------------------------------------------------------------
// tuple<double,double,double,double,double,double,double,double,double> IERS (Matrix eop,double Mjd_UTC,char interp ='n')
//------------------------------------------------------------------------------
/**
 * @brief Operacion IERS.
 *
 * @param [in] eop Vector. 
 * @param [in] Mjd_UTC Tiempo. 
 * @param [in] interp intercepcion. 
 * @return Gestión de datos de tiempo y movimiento polar del IERS.
 */
//------------------------------------------------------------------------------
tuple<double,double,double,double,double,double,double,double,double> IERS (Matrix eop,double Mjd_UTC,char interp ='n');

#endif