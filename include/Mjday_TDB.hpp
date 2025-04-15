//$Header$
//------------------------------------------------------------------------------
// Mjday_TDB
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Mjday_TDB.hpp
 * @brief Archivo cabecera de operacion Mjday_TDB.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _Mjday_TDB_
#define _Mjday_TDB_

#include <cmath>

//------------------------------------------------------------------------------
// double Mjday_TDB(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * @brief Operacion Mjday.
 *
 * @param [in] Mjd_TT Fecha juliana modificada (TT).   
 * @return Fecha juliana modificada (TDB).
 */
//------------------------------------------------------------------------------
double Mjday_TDB(double Mjd_TT);

#endif