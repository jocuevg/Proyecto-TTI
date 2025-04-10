//$Source$
//------------------------------------------------------------------------------
// Cheb3D
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Cheb3D.cpp
 * @brief Programacion de operacion Cheb3D.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#include "..\include\Cheb3D.hpp"

Matrix Cheb3D(double t, double N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz)
{
    if ((t < Ta) || (Tb < t))
        error('ERROR: Time out of range in Cheb3D::Value\n');

    double tau = (2 * t - Ta - Tb) / (Tb - Ta);

    Matrix f1 = zeros(1, 3);
    Matrix f2 = zeros(1, 3);

    for (i = N : -1 : 2)
    {
        Matrix old_f1 = f1;

        Matrix aux = new Matrix(3);
        aux(1) = Cx(i);
        aux(2) = Cy(i);
        aux(3) = Cz(i);

        f1 = f1 * 2 * tau - f2 + aux;
        f2 = old_f1;
    }

    Matrix aux2 = new Matrix(3);
    aux2(1) = Cx(1);
    aux2(2) = Cy(1);
    aux2(3) = Cz(1);

    return f1 * tau - f2 + aux2;
}
