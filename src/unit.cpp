//$Source$
//------------------------------------------------------------------------------
// unit
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file unit.cpp
* @brief Programacion de operacion unit.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\unit.hpp"

Matrix& unit(Matrix& vec){
    double small = 0.000001;
    double magv = norm(transpose(vec));

    Matrix& outvec=zeros(3);

    if ( magv > small ){
        for (int i=1;i<=3;i++){
            outvec(i)= vec(i,1)/magv;
        }
    }
    
    return outvec;
}
