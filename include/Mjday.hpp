//$Header$
//------------------------------------------------------------------------------
// Mjday
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Mjday.hpp
 * @brief Archivo cabecera de operacion Mjday.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _Mjday_
#define _Mjday_

#include <cmath>

//------------------------------------------------------------------------------
// double Mjday(int yr, int mon, int day, int hr, int min, double sec)
//------------------------------------------------------------------------------
/**
 * @brief Operacion Mjday.
 *
 * @param [in] yr Año de la fecha.  
 * @param [in] mon Mes de la fecha.  
 * @param [in] day Día de la fecha.  
 * @param [in] hr Horas de la fecha.  
 * @param [in] min Minutos de la fecha.  
 * @param [in] sec Segundos de la fecha.   
 * @return Fecha juliana modificada.
 */
//------------------------------------------------------------------------------
double Mjday(int yr, int mon, int day, int hr, int min, double sec);

//------------------------------------------------------------------------------
// double Mjday(int yr, int mon, int day)
//------------------------------------------------------------------------------
/**
 * @brief Operacion Mjday.
 *
 * @param [in] yr Año de la fecha.  
 * @param [in] mon Mes de la fecha.  
 * @param [in] day Día de la fecha.    
 * @return Fecha juliana modificada.
 */
//------------------------------------------------------------------------------
double Mjday(int yr, int mon, int day);

#endif