//$Source$
//------------------------------------------------------------------------------
// EKF_GEOS3
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file EKF_GEOS3.cpp
* @brief Initial Orbit Determination using Gauss and Extended Kalman Filter methods.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\Sat_const.hpp"
#include "..\include\Position.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\DEInteg.hpp"
#include "..\include\Accel.hpp"
#include "..\include\LTC.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\gmst.hpp"
#include "..\include\R_z.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\MeasUpdate.hpp"
#include <iostream>
using namespace std;

int main() {

    int nobs = 46;

    eop19620101(21413);
	DE430Coeff(2285,1020);
	GGM03S(181);
	GEOS3(nobs);

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Rad;       // [rad]
    double sigma_el = 0.0139*Rad;       // [rad]

    // Kaena Point station
    double lat = Rad*21.5748;           // [rad]
    double lon = Rad*(-158.2706);       // [rad]
    double alt = 300.20;                // [m]

    Matrix& Rs = transpose(Position(lon, lat, alt));

    double Mjd1 = obs(1,1);
    double Mjd2 = obs(9,1);
    double Mjd3 = obs(18,1);

    Matrix& r2=zeros(3),&v2=zeros(3);
    r2(1)=6221397.62857869;r2(2)=2867713.77965738;r2(3)=3006155.98509949;
    v2(1)=4645.04725161806;v2(2)=-2752.21591588204;v2(3)=-7507.99940987031;
    // [r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
    // [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
    //                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

    Matrix& Y0_apr = transpose(union_vector(r2,v2));

    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n      = 20;
    AuxParam.m      = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;

    int n_eqn  = 6;

    double t_aux=0.0;

    Matrix& Y = DEInteg(Accel,t_aux,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);

    Matrix& P = zeros(6,6);
    
    for (int i=1;i<=3;i++){
        P(i,i)=1e8;
    }
    for (int i=4;i<=6;i++){
        P(i,i)=1e3;
    }

    Matrix& LT = LTC(lon,lat);

    Matrix& yPhi = zeros(42,1);
    Matrix& Phi  = zeros(6,6);

    // Measurement loop
    double t = 0.0;

    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,Mjd_TT,Mjd_UT1,t_old,theta,Azim, Elev;

    for (int i=1;i<=nobs;i++){ 
        
        // Previous step
        t_old = t;
        Matrix& Y_old = Y;
        
        // Time increment and propagation
        Mjd_UTC = obs(i,1);                       // Modified Julian Date
        t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
        
        tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(Mjd_UTC,'l');
        tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC/86400;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;
            
        for (int ii=1;ii<=6;ii++){
            yPhi(ii,1) = Y_old(ii,1);
            for (int j=1;j<=6;j++){  
                if (ii==j){ 
                    yPhi(6*j+ii,1) = 1; 
                }else{
                    yPhi(6*j+ii,1) = 0;
                }
            }
        }

        cout<<i<<":\n\n\n\n";
        yPhi = DEInteg(VarEqn,t_aux,t-t_old,1e-13,1e-6,42,yPhi);
        cout<<yPhi<<"\n";

        // Extract state transition matrices
        for (int j=1;j<=6;j++){
            assign_column(Phi,extract_vector(transpose(yPhi),6*j+1,6*j+6),j);
        }

        cout<<"Fallo:\n\n";
        Y = DEInteg(Accel,t_aux,t-t_old,1e-13,1e-6,6,Y_old);
        cout<<Y<<"\n";
        
        // Topocentric coordinates
        theta = gmst(Mjd_UT1);                    // Earth rotation
        Matrix& U = R_z(theta);
        Matrix& r = transpose(extract_vector(transpose(Y),1,3));
        Matrix& s = LT*(U*r-Rs);                          // Topocentric position [m]
        
        // Time update
        P = TimeUpdate(P, Phi);
            
        // Azimuth and partials
        Matrix &dAds=zeros(3),&dEds=zeros(3);
        tie(Azim, Elev, dAds, dEds) = AzElPa(s);     // Azimuth, Elevation
        Matrix& dAdY = union_vector(dAds*LT*U,zeros(1,3));
        
        // Measurement update
        Matrix &K=zeros(6,1);

        tie(K, Y, P) = MeasUpdate ( Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );
        
        // Elevation and partials
        r = transpose(extract_vector(transpose(Y),1,3));
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        tie(Azim, Elev, dAds, dEds) = AzElPa(s);     // Azimuth, Elevation
        Matrix& dEdY = union_vector(dEds*LT*U,zeros(1,3));
        
        // Measurement update
        tie(K, Y, P) = MeasUpdate ( Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );
        
        // Range and partials
        r = transpose(extract_vector(transpose(Y),1,3));
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        double Dist = norm(transpose(s)); Matrix& dDds = transpose(s/Dist);         // Range
        Matrix& dDdY = union_vector(dDds*LT*U,zeros(1,3));
        
        // Measurement update
        tie(K, Y, P) = MeasUpdate ( Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );

    }

    cout<<"AA\n";

    tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(obs(46,1),'l');
    tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix& Y0 = DEInteg(Accel,t_aux,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

    Matrix& Y_true = zeros(6); 
    Y_true(1)= 5753.173e3; 
    Y_true(2)= 2673.361e3;
    Y_true(3)= 3440.304e3;
    Y_true(4)= 4.324207e3;
    Y_true(5)= -1.924299e3;
    Y_true(6)= -5.728216e3;

    printf("\nError of Position Estimation\n");
    printf("dX%10.1f [m]\n",Y0(1,1)-Y_true(1));
    printf("dY%10.1f [m]\n",Y0(2,1)-Y_true(2));
    printf("dZ%10.1f [m]\n",Y0(3,1)-Y_true(3));
    printf("\nError of Velocity Estimation\n");
    printf("dVx%8.1f [m/s]\n",Y0(4,1)-Y_true(4));
    printf("dVy%8.1f [m/s]\n",Y0(5,1)-Y_true(5));
    printf("dVz%8.1f [m/s]\n",Y0(6,1)-Y_true(6));

    return 0;
}