//$Source$
//------------------------------------------------------------------------------
// angl
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file angl.cpp
* @brief Programacion de operacion angl.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\angl.hpp"
#include "..\include\sign_.hpp"

double angl(Matrix& vec1,Matrix& vec2){

    double small     = 0.00000001;
    double undefined = 999999.1;

    double magv1 = norm(transpose(vec1));
    double magv2 = norm(transpose(vec2));

    double theta;

    if (magv1*magv2 > small^2){
        double temp= dot(transpose(vec1),transpose(vec2)) / (magv1*magv2);
        if (abs( temp ) > 1.0){
            int sign;
            if (temp > 0) sign= 1;
            else if (temp < 0) sign= -1;
            else sign= 0;
            temp= sign * 1.0;
        }
        theta= acos( temp );
    }else{
        theta= undefined;
    }

    return theta;
}
