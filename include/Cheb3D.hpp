//$Header$
//------------------------------------------------------------------------------
// Cheb3D
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Cheb3D.hpp
 * @brief Archivo cabecera de operacion Cheb3D.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _Cheb3D_
#define _Cheb3D_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& Cheb3D(double t,double N,double Ta,double Tb,Matrix & Cx,Matrix& Cy,Matrix& Cz)
//------------------------------------------------------------------------------
/**
 * @brief Operacion Cheb3D.
 *
 * @param [in] t Tiempo.
 * @param [in] N Numero de coeficiente.
 * @param [in] Ta Inicio del intervalo.
 * @param [in] Tb Final del intervale.
 * @param [in] Cx Coeficientes del polimonio de Chebyshev (coordenada x) fila.
 * @param [in] Cy Coeficientes del polimonio de Chebyshev (coordenada y) fila.
 * @param [in] Cz Coeficientes del polimonio de Chebyshev (coordenada z) fila.
 * @return Vector fila, aproximación de Chebyshev de vectores de 3 dimensiones.
 */
//------------------------------------------------------------------------------
 Matrix& Cheb3D(double t,double N,double Ta,double Tb,Matrix& Cx,Matrix& Cy,Matrix& Cz);

#endif