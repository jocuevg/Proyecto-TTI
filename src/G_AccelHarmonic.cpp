//$Source$
//------------------------------------------------------------------------------
// G_AccelHarmonic
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file G_AccelHarmonic.cpp
 * @brief Programacion de operacion G_AccelHarmonic.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\AccelHarmonic.hpp"

Matrix& G_AccelHarmonic(Matrix &r, Matrix &U, int n_max, int m_max)
{
    double d = 1.0; // Position increment [m]

    Matrix &G = zeros(3, 3);
    Matrix &dr = zeros(3, 1);

    // Gradient
    for (int i = 1; i <= 3; i++)
    {
        // Set offset in i - th component of the position vector
        dr=zeros(3,1);
        dr(i,1) = d;
        // Acceleration difference
        Matrix& da = AccelHarmonic(r + dr / 2, U, n_max, m_max) - AccelHarmonic(r - dr / 2, U, n_max, m_max);
        // Derivative with respect to i - th axis
        assign_column(G,transpose(da/d),i);
    }

    return G;
}
