//$Source$
//------------------------------------------------------------------------------
// AzElPa
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file AzElPa.cpp
 * @brief Programacion de operacion AzElPa.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#include "..\include\AzElPa.hpp"
#include "..\include\Sat_const.hpp"
using namespace std;

tuple<double, double, Matrix, Matrix> AzElPa(Matrix s)
{
    double pi2 = 2.0 * pi;

    double rho = sqrt(s(1) * s(1) + s(2) * s(2));

    // Angles
    double Az = atan2(s(1), s(2));

    if (Az < 0.0)
        Az = Az + pi2;

    double El = atan(s(3) / rho);

    // Partials
    Matrix dAds(3);
    dAds(1) = s(2) / (rho * rho);
    dAds(2) = -s(1) / (rho * rho);
    dAds(3) = 0;

    Matrix dEds(3);
    dEds(1) = -s(1) * s(3) / rho;
    dEds(2) = -s(2) * s(3) / rho;
    dEds(3) = rho;
    dEds = dEds / dot(s, s);

    return make_tuple(Az, El, dAds, dEds);
}
