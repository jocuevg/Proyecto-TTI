//$Header$
//------------------------------------------------------------------------------
// Position
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Position.hpp
 * @brief Archivo cabecera de operacion Position.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _Position_
#define _Position_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& Position(double lon,double lat,double h)
//------------------------------------------------------------------------------
/**
 * @brief Operacion Position.
 *
 * @param [in] lon Longitud en radianes.   
 * @param [in] lan Latitud en radianes.   
 * @param [in] h Altura en metros.   
 * @return Vector fila posicion.
 */
//------------------------------------------------------------------------------
Matrix& Position(double lon,double lat,double h);

#endif