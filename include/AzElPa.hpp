//$Header$
//------------------------------------------------------------------------------
// AzElPa
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file AzElPa.hpp
 * @brief Archivo cabecera de operacion AzElPa.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _AzElPa_
#define _AzElPa_

#include <cmath>
#include <tuple>
#include "matrix.hpp"
using namespace std;

//------------------------------------------------------------------------------
// tuple<double,double,Matrix,Matrix> AzElPa (Matrix s)
//------------------------------------------------------------------------------
/**
 * @brief Operacion AzElPa.
 *
 * @param [in] s Coordenadas de la tangente local topocéntrica. 
 * @return Acimut [rad], Elevación [rad], Parciales de acimut con respecto a s , y Parciales de elevación con respecto a s.
 */
//------------------------------------------------------------------------------
tuple<double,double,Matrix,Matrix> AzElPa (Matrix s);

#endif