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

Matrix& Cheb3D(double t, double N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz)
{
    if ((t < Ta) || (Tb < t)){
        cout << "ERROR: Time out of range in Cheb3D::Value\n";
		exit(EXIT_FAILURE);
    }

    double tau = (2 * t - Ta - Tb) / (Tb - Ta);

    Matrix& f1 = zeros(3);
    Matrix& f2 = zeros(3);

    for (int i=N;i>=2;i--)
    {
        Matrix& old_f1=zeros(3);
        old_f1=f1;

        Matrix& aux=zeros(3);
        aux(1) = Cx(i);
        aux(2) = Cy(i);
        aux(3) = Cz(i);

        f1 = f1 * 2 * tau - f2 + aux;

        f2 = old_f1;
    }

    Matrix& aux2=zeros(3);
    aux2(1) = Cx(1);
    aux2(2) = Cy(1);
    aux2(3) = Cz(1);

    return f1 * tau - f2 + aux2;
}
