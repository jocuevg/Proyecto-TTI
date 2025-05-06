//$Source$
//------------------------------------------------------------------------------
// MeasUpdate
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file MeasUpdate.cpp
* @brief Programacion de operacion MeasUpdate.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\MeasUpdate.hpp"

tuple<Matrix&,Matrix&,Matrix&> MeasUpdate(Matrix& x,double z,double g,double s,Matrix& G,Matrix& P,int n){

    int m = 1;
    Matrix& Inv_W = zeros(m,m);

    for (int i=1;i<=m;i++){
        Inv_W(i,i) = s*s;    // Inverse weight (measurement covariance)
    }

    // Kalman gain
    Matrix& K = P*transpose(G)*inv(Inv_W+G*P*transpose(G))(1,1);

    // State update
    x = x + K*(z-g);

    // Covariance update
    P = (eye(n)-K*G)*P;

    return tie(K,x,P);

}
