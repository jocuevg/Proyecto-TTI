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

tuple<double, double, Matrix&, Matrix&> AzElPa(Matrix& s)
{
    double pi2 = 2.0 * pi;

    Matrix& st=transpose(s);

    double rho = sqrt(st(1) * st(1) + st(2) * st(2));

    // Angles
    double Az = atan2(st(1), st(2));

    if (Az < 0.0)
        Az = Az + pi2;

    double El = atan(st(3) / rho);

    // Partials
    Matrix &dAds = zeros(3);
    dAds(1) = st(2) / (rho * rho);
    dAds(2) = -st(1) / (rho * rho);
    dAds(3) = 0;

    Matrix &dEds = zeros(3);
    dEds(1) = -st(1) * st(3) / rho;
    dEds(2) = -st(2) * st(3) / rho;
    dEds(3) = rho;
    dEds = dEds / dot(st, st);

    return tie(Az, El, dAds, dEds);
}
