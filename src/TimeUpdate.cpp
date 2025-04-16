//$Source$
//------------------------------------------------------------------------------
// TimeUpdate
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file TimeUpdate.cpp
 * @brief Programacion de operacion TimeUpdate.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#include "..\include\TimeUpdate.hpp"
#include "..\include\Sat_const.hpp"
using namespace std;

Matrix TimeUpdate(Matrix P, Matrix Phi, double Qdt)
{
  return Phi * P * transpose(Phi) + Qdt;
}
