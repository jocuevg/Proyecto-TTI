//$Header$
//------------------------------------------------------------------------------
// VarEqn
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file VarEqn.hpp
 * @brief Archivo cabecera de operacion VarEqn.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _VarEqn_
#define _VarEqn_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& VarEqn(double x,Matrix& yPhi)
//------------------------------------------------------------------------------
/**
 * @brief Operacion VarEqn.
 *
 * @param [in] x Tiempo desde epoch en segundos.
 * @param [in] Y vector columna (6+36)-dim que comprende el vector de estado (y) 
 *      y la matriz de transición de estado (Phi) en orden de almacenamiento por columnas.
 * @return Derivada de yPhi (vector fila por conveniencia).
 */
//------------------------------------------------------------------------------
Matrix& VarEqn(double x,Matrix& yPhi);

#endif