//$Header$
//------------------------------------------------------------------------------
// TimeUpdate
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file TimeUpdate.hpp
 * @brief Archivo cabecera de operacion TimeUpdate.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _TimeUpdate_
#define _TimeUpdate_

#include <cmath>
#include "matrix.hpp"
using namespace std;

//------------------------------------------------------------------------------
// Matrix& TimeUpdate(Matrix& P, Matrix& Phi, double Qdt = 0.0)
//------------------------------------------------------------------------------
/**
 * @brief Operacion TimeUpdate.
 *
 * @param [in] P Hora original.
 * @param [in] Phi Nueva.
 * @param [in] Qdt Cuadrante.
 * @return Hora adaptada.
 */
//------------------------------------------------------------------------------
Matrix& TimeUpdate(Matrix& P, Matrix& Phi, double Qdt = 0.0);

#endif