//$Header$
//------------------------------------------------------------------------------
// DEInteg
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file DEInteg.hpp
 * @brief Archivo cabecera de operacion DEInteg.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _DEInteg_
#define _DEInteg_

#include <cmath>
#include "matrix.hpp"

//------------------------------------------------------------------------------
// Matrix& DEInteg(Matrix& f (double t, Matrix& y), double t, double tout, double relerr, double abserr,int n_eqn, Matrix& y)
//------------------------------------------------------------------------------
/**
 * @brief Operacion DEInteg.
 *
 * @param [in] func Funcion de integracion.
 * @param [in] t Inicio 
 * @param [in] tout Fin
 * @param [in] relerr Error relativo 
 * @param [in] abserr Error absoluto
 * @param [in] n_eqn Numero de ecuaciones 
 * @param [in] y Vector columna de coeficientes
 * @return Integral de y (vector columna).
 */
//------------------------------------------------------------------------------
Matrix& DEInteg(Matrix& func (double t, Matrix& y), double t, double tout, double relerr, double abserr,int n_eqn, Matrix& y);

typedef struct {
    int DE_INIT,DE_DONE,DE_BADACC,DE_NUMSTEPS,DE_STIFF,DE_INVPARAM;
} STATE;

//------------------------------------------------------------------------------
// STATE DE_STATE
//------------------------------------------------------------------------------
/**
 * @brief Variable DE_STATE.
 */
//------------------------------------------------------------------------------
extern STATE DE_STATE;

#endif