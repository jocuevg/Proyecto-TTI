//$Header$
//------------------------------------------------------------------------------
// PoleMatrix
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file PoleMatrix.hpp
 * @brief Archivo cabecera de operacion PoleMatrix.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _PoleMatrix_
#define _PoleMatrix_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& PoleMatrix(double xp, double yp)
//------------------------------------------------------------------------------
/**
 * @brief Operacion PoleMatrix.
 *
 * @param [in] xp Coordenada x del polo.
 * @param [in] yp Coordenada y del polo.
 * @return Matriz del polo.
 */
//------------------------------------------------------------------------------
Matrix& PoleMatrix(double xp, double yp);

#endif