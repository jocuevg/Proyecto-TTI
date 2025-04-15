//$Source$
//------------------------------------------------------------------------------
// EccAnom
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file EccAnom.cpp
 * @brief Programacion de operacion EccAnom.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#include "..\include\EccAnom.hpp"
#include "..\include\Sat_const.hpp"

double EccAnom(double M, double e)
{
    int maxit = 15;
    int i = 1;
    double E;

    // Starting value
    M = fmod(M, 2.0 * pi);

    if (e < 0.8)
        E = M;
    else
        E = pi;

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    // Iteration
    while (fabs(f) > 1e2 * __DBL_EPSILON__)
    {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i = i + 1;
        if (i == maxit)
        {
            cout << " convergence problems in EccAnom";
            exit(EXIT_FAILURE);
        }
    }
    return E;
}
