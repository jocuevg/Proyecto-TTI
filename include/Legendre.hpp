//$Header$
//------------------------------------------------------------------------------
// Legendre
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Legendre.hpp
 * @brief Archivo cabecera de operacion Legendre.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _Legendre_
#define _Legendre_

#include <cmath>
#include <tuple>
#include "matrix.hpp"
using namespace std;

//------------------------------------------------------------------------------
// tuple<Matrix&,Matrix&> Legendre (int n,int m,double fi)
//------------------------------------------------------------------------------
/**
 * @brief Operacion Legendre.
 *
 * @param [in] n Fila-1. 
 * @param [in] m Columnas-1. 
 * @param [in] fi Radianes. 
 * @return Polinomio de Legendre y derivado.
 */
//------------------------------------------------------------------------------
tuple<Matrix&,Matrix&> Legendre (int n,int m,double fi);

#endif