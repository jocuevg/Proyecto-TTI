//$Source$
//------------------------------------------------------------------------------
// tests
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file tests.cpp
* @brief Tests de funciones.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\global.hpp"
#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\Sat_const.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\Position.hpp"
#include "..\include\sign_.hpp"
#include "..\include\timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\IERS.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\LTC.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\gmst.hpp"
#include "..\include\gast.hpp"
#include "..\include\MeasUpdate.hpp"
#include <cstdio>
#include <cmath>
#include <iomanip>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

bool compareTuplesAzElPa(const tuple<double, double, Matrix, Matrix> t1, const tuple<double, double, Matrix, Matrix> t2, double p) {
	if(fabs(get<0>(t1)-get<0>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<0>(t1),get<0>(t2));
		return 0;
	}
	if(fabs(get<1>(t1)-get<1>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<1>(t1),get<1>(t2));
		return 0;
	}

    return m_equals(get<2>(t1), get<2>(t2),p) && m_equals(get<3>(t1), get<3>(t2),p);
}

bool compareTuplesIERS(tuple<double, double, double, double, double, double, double, double, double> t1, 
		tuple<double, double, double, double, double, double, double, double, double> t2, 
		double p) {

	if(fabs(get<0>(t1)-get<0>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<0>(t1),get<0>(t2));
		return 0;
	}
	if(fabs(get<1>(t1)-get<1>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<1>(t1),get<1>(t2));
		return 0;
	}
	if(fabs(get<2>(t1)-get<2>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<2>(t1),get<2>(t2));
		return 0;
	}
	if(fabs(get<3>(t1)-get<3>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<3>(t1),get<3>(t2));
		return 0;
	}
	if(fabs(get<4>(t1)-get<4>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<4>(t1),get<4>(t2));
		return 0;
	}
	if(fabs(get<5>(t1)-get<5>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<5>(t1),get<5>(t2));
		return 0;
	}
	if(fabs(get<6>(t1)-get<6>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<6>(t1),get<6>(t2));
		return 0;
	}
	if(fabs(get<7>(t1)-get<7>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<7>(t1),get<7>(t2));
		return 0;
	}
	if(fabs(get<8>(t1)-get<8>(t2)) > p) {
		printf("%2.20lf %2.20lf\n",get<8>(t1),get<8>(t2));
		return 0;
	}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_zeros_02() {
    int f = 3;
	
	Matrix A(f);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0;

	Matrix B = zeros(3);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_mul_01() {
    int f = 3;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; 
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; 
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; 
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; 
	
	Matrix C(f, c);
	C(1,1) = 14; C(1,2) =  -28; C(1,3) = 2; 
	C(2,1) = -5; C(2,2) = 2; C(2,3) = -1; 
	C(3,1) = 7; C(3,2) = -2; C(3,3) = 1; 
	
	Matrix R = A * B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_div_01() {
    int f = 3;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; 
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; 
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; 
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; 
	
	Matrix C(f, c);
	C(1,1) = -28; C(1,2) = 8; C(1,3) = -6; 
	C(2,1) = 0.5; C(2,2) = 0; C(2,3) = 0.3333333333333333; 
	C(3,1) = 0; C(3,2) = 0; C(3,3) = -0.3333333333333333; 
	
	Matrix R = A / B;

    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_eq_01() {
    int f = 3;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; 
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; 
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; 
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; 
	
	Matrix C(f, c);
	C(1,1) = 0; C(1,2) =  2; C(1,3) = 8;
	C(2,1) = 1; C(2,2) = -1; C(2,3) = 0;
	C(3,1) = 0; C(3,2) =  1; C(3,3) = 0; 
	
	B = A ;
    
    _assert(m_equals(C, B, 1e-10));
    
    return 0;
}

int m_eye_01() {
    int f = 3;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0; 
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1; 
	
	Matrix B = eye(3);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_trans_01() {
    int f = 3;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 1; 
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1; 
	
	Matrix B(f, c);
	B(1,1) = 1; B(1,2) = 0; B(1,3) = 0; 
	B(2,1) = 0; B(2,2) = 1; B(2,3) = 0; 
	B(3,1) = 1; B(3,2) = 0; B(3,3) = 1;

	Matrix R=transpose(A);
    
    _assert(m_equals(R, B, 1e-10));
    
    return 0;
}

int m_inv_01() {
    int f = 3;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0;  
	
	Matrix B(f, c);
	B(1,1) = 0; B(1,2) = 1; B(1,3) = 1; 
	B(2,1) = 0; B(2,2) = 0; B(2,3) = 1; 
	B(3,1) = 0.125; B(3,2) = 0; B(3,3) = -0.25;

	Matrix R=inv(A);
    
    _assert(m_equals(R, B, 1e-10));
    
    return 0;
}

int m_sum_02() {
    int f = 3;
    int c = 4;
	double x = 1.5;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix C(f, c);
	C(1,1) = 1.5; C(1,2) =  3.5; C(1,3) = 9.5; C(1,4) = 1.5;
	C(2,1) = 2.5; C(2,2) = 0.5; C(2,3) = 1.5; C(2,4) = 1.5;
	C(3,1) = 1.5; C(3,2) = 2.5; C(3,3) = 1.5; C(3,4) = 6.5;
	
	Matrix R = A + x;

    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_02() {
    int f = 3;
    int c = 4;
	double x = 1.5;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix C(f, c);
	C(1,1) = -1.5; C(1,2) =  0.5; C(1,3) = 6.5; C(1,4) = -1.5;
	C(2,1) = -0.5; C(2,2) = -2.5; C(2,3) = -1.5; C(2,4) = -1.5;
	C(3,1) = -1.5; C(3,2) = -0.5; C(3,3) = -1.5; C(3,4) = 3.5;
	
	Matrix R = A - x;

    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_mul_02() {
    int f = 3;
    int c = 4;
	double x = 1.5;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix C(f, c);
	C(1,1) = 0; C(1,2) =  3; C(1,3) = 12; C(1,4) = 0;
	C(2,1) = 1.5; C(2,2) = -1.5; C(2,3) = 0; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 1.5; C(3,3) = 0; C(3,4) = 7.5;
	
	Matrix R = A * x;

    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_div_02() {
    int f = 3;
    int c = 4;
	double x = 2.5;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix C(f, c);
	C(1,1) = 0; C(1,2) =  0.8; C(1,3) = 3.2; C(1,4) = 0;
	C(2,1) = 0.4; C(2,2) = -0.4; C(2,3) = 0; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 0.4; C(3,3) = 0; C(3,4) = 2;
	
	Matrix R = A / x;

    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_norm_01() {
    int f = 3;

	Matrix A(f);
	A(1,1) = 0; A(1,2) =  3; A(1,3) = 4;


    _assert(norm(A)==5);
    
    return 0;
}

int m_dot_01() {
    int f = 3;

	Matrix A(f);
	A(1,1) = 0; A(1,2) =  3; A(1,3) = 4;

	Matrix B(f);
	B(1,1) = 0; B(1,2) =  2; B(1,3) = 1;


    _assert(dot(A,B)==10);
    
    return 0;
}

int m_cross_01() {
    int f = 3;

	Matrix A(f);
	A(1,1) = 0; A(1,2) =  3; A(1,3) = 4;

	Matrix B(f);
	B(1,1) = 0; B(1,2) =  2; B(1,3) = 1;

	Matrix R(f);
	R(1,1) = -5; R(1,2) =  0; R(1,3) = 0;

	Matrix C=cross(A,B);

	_assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_extractVector_01() {
    int c = 4;
	
	Matrix A(c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;

	Matrix B(2);
	B(1,1) = 2; B(1,2) =  8; 

	Matrix R=extract_vector(A,2,3);

	_assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_extractRow_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;

	Matrix B(c);
	B(1,1) = 0; B(1,2) =  2; B(1,3) = 8; B(1,4) = 0;

	Matrix R=extract_row(A,1);

	_assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_extractColumn_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;

	Matrix B(f);
	B(1,1) = 2; B(1,2) =  -1; B(1,3) = 1;

	Matrix R=extract_column(A,2);

	_assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_unionVector_01() {
    int c = 4;
	
	Matrix A(c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;

	Matrix B(c);
	B(1,1) = 0; B(1,2) =  2; B(1,3) = 8; B(1,4) = 0;

	Matrix C(2*c);
	C(1,1) = 0; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0; C(1,5) = 0; C(1,6) =  2; C(1,7) = 8; C(1,8) = 0;

	Matrix R=union_vector(A,B);

	_assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_assignRow_01() {
	int f = 3;
    int c = 3;

    Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0;

	Matrix B(c);
	B(1,1) = 99; B(1,2) =  100; B(1,3) = 101;

	Matrix C(f,c);
	C(1,1) = 99; C(1,2) =  100; C(1,3) = 101;
	C(2,1) = 1; C(2,2) = -1; C(2,3) = 0;
	C(3,1) = 0; C(3,2) =  1; C(3,3) = 0;

	Matrix R=assign_row(A,B,1);

	_assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_assignColumn_01() {
    int f = 3;
    int c = 3;

    Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0;

	Matrix B(c);
	B(1,1) = 99; B(1,2) =  100; B(1,3) = 101;

	Matrix C(f,c);
	C(1,1) = 0; C(1,2) =  99; C(1,3) = 8;
	C(2,1) = 1; C(2,2) = 100; C(2,3) = 0;
	C(3,1) = 0; C(3,2) = 101; C(3,3) = 0;

	Matrix R=assign_column(A,B,2);

	_assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int R_x_01() {
    int f = 3;
    int c = 3;

    Matrix A(f, c);
	A(1,1) = 1; A(1,2) =  0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = -0.989992496600445; A(2,3) = 0.141120008059867;
	A(3,1) = 0; A(3,2) =  -0.141120008059867; A(3,3) = -0.989992496600445;

	Matrix R=R_x(3);

	_assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int R_y_01() {
    int f = 3;
    int c = 3;

    Matrix A(f, c);
	A(1,1) = -0.989992496600445; A(1,2) =  0; A(1,3) = -0.141120008059867;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0.141120008059867; A(3,2) =  0; A(3,3) = -0.989992496600445;

	Matrix R=R_y(3);

	_assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int R_z_01() {
    int f = 3;
    int c = 3;

    Matrix A(f, c);
	A(1,1) = -0.989992496600445; A(1,2) =  0.141120008059867; A(1,3) = 0;
	A(2,1) = -0.141120008059867; A(2,2) = -0.989992496600445; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  0; A(3,3) = 1;

	Matrix R=R_z(3);

	_assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int AccelPointMass_01() {
    int f = 3;
    int c = 3;
	
	Matrix A(f,1);
	A(1,1) = 1; A(2,1) =  2; A(3,1) = 3; 
	
	Matrix B(f,1);
	B(1,1) = 4; B(2,1) =  5; B(3,1) = 6;
	
	Matrix C(f,1);
	C(1,1) = 0.0463899400720928; C(2,1) =  0.0419499176126264; C(3,1) = 0.0375098951531599; 
	
	Matrix R = AccelPointMass(A,B,3);
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int Cheb3D_01() {
    int f = 3;
	
	Matrix A(f);
	A(1,1) = 1; A(1,2) =  2; A(1,3) = 3; 
	
	Matrix B(f);
	B(1,1) = 4; B(1,2) =  5; B(1,3) = 6;
	
	Matrix C(f);
	C(1,1) = 7; C(1,2) =  8; C(1,3) = 9; 

	Matrix D(f);
	D(1,1) = -2; D(1,2) =  -2; D(1,3) = -2; 
	
	Matrix R = Cheb3D(5,3,0,10,A,B,C);
    
    _assert(m_equals(D, R, 1e-10));
    
    return 0;
}

int EccAnom_01() {

    _assert(fabs(2.3542427582227807 - EccAnom(2,0.5)) < 1e-10);
    
    return 0;
}

int Frac_01() {

    _assert(Frac(2.5)==0.5);
    
    return 0;
}

int MeanObliquity_01() {

    _assert(fabs(0.4094130391 - MeanObliquity(5)) < 1e-10);
    
    return 0;
}

int Mjday_01(){

	_assert(fabs(52347.4989583333 - Mjday(2002,3,14,11,58,30)) < 1e-10);
    
    return 0;
}

int Mjday_02(){

	_assert(52347 == Mjday(2002,3,14));
    
    return 0;
}

int MjdayTDB_01(){

	_assert(fabs(52347.4989583509 - Mjday_TDB(52347.4989583333)) < 1e-10);
    
    return 0;
}

int Position_01(){

	Matrix A(3);
	A(1)=219206.013888262; 
	A(2)=-5500434.57124668; 
	A(3)=-3210856.41301074;

	Matrix R= Position(99,100,101);

	_assert(m_equals(A, R, 1e-8));
    
    return 0;
}

int sign_01(){

	_assert(-8.5==sign_(8.5,-1));
    
    return 0;
}

int timediff_01(){

	tuple<double,double,double,double,double> A = timediff(5,10);
	
	tuple<double,double,double,double,double> B = make_tuple(-5,9,14,42.184,9);
	
	_assert(A==B);

	return 0;
}

int AzElPa_01(){

	Matrix A(3,1);
	A(1,1)=1;A(2,1)=1;A(3,1)=1;
	tuple<double,double,Matrix&,Matrix&> R = AzElPa(A);
	
	Matrix& B=zeros(3);
	B(1)=0.5;B(2)=-0.5;B(3)=0;
	Matrix& C=zeros(3);
	C(1)=-0.23570226039;C(2)=-0.23570226039;C(3)=0.47140452079;

	double x=0.78539816339;
	double y=0.61547970867;
	tuple<double,double,Matrix&,Matrix&> D = tie(x,y,B,C);
	
	_assert(compareTuplesAzElPa(R,D,1e-10));

	return 0;
}

int IERS_01(){

	tuple<double, double, double, double, double, double, double, double, double> R = IERS(47777,'l');

	tuple<double, double, double, double, double, double, double, double, double> D = 
		make_tuple(1.11028150738257e-06,7.70679220038963e-07,-0.4566596,0.0006642,-6.12222716505122e-08,-1.22560898584491e-08,1.93925472443814e-09,3.10280755910103e-10,24);

	_assert(compareTuplesIERS(R,D,1e-10));

	return 0;
}

int IERS_02(){

	tuple<double, double, double, double, double, double, double, double, double> R = IERS(47777);

	tuple<double, double, double, double, double, double, double, double, double> D = 
		make_tuple(1.11028150738257e-06,7.70679220038963e-07,-0.4566596,0.0006642,-6.12222716505122e-08,-1.22560898584491e-08,1.93925472443814e-09,3.10280755910103e-10,24);

	_assert(compareTuplesIERS(R,D,1e-10));

	return 0;
}

int Legendre_01(){

	Matrix A(3,3);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0; 
	A(2,1) = 0; A(2,2) = -1.732050807568877; A(2,3) = 0;
	A(3,1) = -1.118033988749895; A(3,2) = 0; A(3,3) = 1.936491673103709; 

	Matrix B(3,3);
	B(1,1) = 0; B(1,2) =  0; B(1,3) = 0; 
	B(2,1) = -1.732050807568877; B(2,2) = 0; B(2,3) = 0; 
	B(3,1) = 0; B(3,2) = 3.872983346207417; B(3,3) = 0;

	tuple<Matrix,Matrix> R = Legendre(2,2,pi);

	_assert(m_equals(A, get<0>(R), 1e-10) && m_equals(B, get<1>(R), 1e-10));

	return 0;
}

int NutAngles_01(){

	tuple<double,double> R = NutAngles(2025);

	_assert((fabs(5.81895359450653e-05 - get<0>(R)) < 1e-10) && (fabs(-3.27016318829544e-05 - get<1>(R)) < 1e-10));
    
    return 0;
}

int TimeUpdate_01() {
    int f = 3;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; 
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; 
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; 
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; 
	
	Matrix C(f, c);
	C(1,1) = 1.5; C(1,2) =  9.5; C(1,3) = -10.5; 
	C(2,1) = -2.5; C(2,2) = 9.5; C(2,3) = -49.5; 
	C(3,1) = -4.5; C(3,2) = -25.5; C(3,3) = -7.5; 
	
	Matrix R = TimeUpdate(A,B,1.5);
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int TimeUpdate_02() {
    int f = 3;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; 
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; 
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; 
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; 
	
	Matrix C(f, c);
	C(1,1) = 0; C(1,2) =  8; C(1,3) = -12; 
	C(2,1) = -4; C(2,2) = 8; C(2,3) = -51; 
	C(3,1) = -6; C(3,2) = -27; C(3,3) = -9; 
	
	Matrix R = TimeUpdate(A,B);
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int JPL_Eph_DE430_01(){
	auto [A,B,C,D,E,F,G,H,I,J,K] = JPL_Eph_DE430(52000);

	Matrix A1(3);
	A1(1)=177822310689.899;A1(2)=-22536323956.0336;A1(3)=-17986999536.9575;
	Matrix B1(3);
	B1(1)=41424286859.9559;B1(2)=3305649882.89728;B1(3)=7701979070.65503;
	Matrix C1(3);
	C1(1)=-147162345749.885;C1(2)=-27647512945.7394;C1(3)=-11963874105.9553;
	Matrix D1(3);
	D1(1)=-21913616871.6504;D1(2)=-122540647338.219;D1(3)=-52326196131.0471;
	Matrix E1(3);
	E1(1)=316858333084.401;E1(2)=708321272026.495;E1(3)=299592115258.39;
	Matrix F1(3);
	F1(1)=778007263812.748;F1(2)=1150810008379.17;F1(3)=448715670551.606;
	Matrix G1(3);
	G1(1)=2479329392309.5;G1(2)=-1669320965559.57;G1(3)=-764251929876.677;
	Matrix H1(3);
	H1(1)=2834369256371.92;H1(2)=-3293467008979.32;H1(3)=-1414288770524.34;
	Matrix I1(3);
	I1(1)=-1122924880805.08;I1(2)=-4230959964564.62;I1(3)=-934343764961.337;
	Matrix J1(3);
	J1(1)=-35792541.6234018;J1(2)=342500830.423682;J1(3)= 145221181.58225;
	Matrix K1(3);
	K1(1)=146580480509.696;K1(2)=26921924222.0793;K1(3)=11671830692.0667;

	_assert(m_equals(A,transpose(A1), 1e-1));
	_assert(m_equals(B,transpose(B1), 1e-1));
	_assert(m_equals(C,transpose(C1), 1e-1));
	_assert(m_equals(D,transpose(D1), 1e-1));
	_assert(m_equals(E,transpose(E1), 1e-1));
	_assert(m_equals(F,transpose(F1), 1e-1));
	_assert(m_equals(G,transpose(G1), 1e-1));
	_assert(m_equals(H,transpose(H1), 1e-1));
	_assert(m_equals(I,transpose(I1), 1e-1));
	_assert(m_equals(J,transpose(J1), 1e-1));
	_assert(m_equals(K,transpose(K1), 1e-1));
	
	return 0;
}

int AccelHarmonic_01(){
	Matrix r(3,1);
	r(1,1)=1;r(2,1)=1;r(3,1)=1;
	Matrix E(3,3);
	E(1,1)=1;E(1,2)=2;E(1,3)=3;
	E(2,1)=4;E(2,2)=5;E(2,3)=6;
	E(3,1)=7;E(3,2)=8;E(3,3)=9;
	Matrix A=AccelHarmonic(r,E,3,3);

	Matrix R(3,1);
	R(1,1)=-1.21291106006315e+23;R(2,1)=-1.4722347029891e+23;R(3,1)=-1.73155834591505e+23;

	_assert(m_equals(A,R, 1e9));
	
	return 0;
}

int EqnEquinox_01(){
	_assert(fabs(EqnEquinox(2025)-5.33807310610183e-05)<1e-10);

	return 0;
}

int LTC_01(){
	Matrix E(3,3);
	E(1,1)=0;E(1,2)=1;E(1,3)=0;
	E(2,1)=0.54402111088937;E(2,2)=0;E(2,3)=-0.839071529076452;
	E(3,1)=-0.839071529076452;E(3,2)=0;E(3,3)=-0.54402111088937;

	_assert(m_equals(E,LTC(0,10), 1e-10));

	return 0;
}

int NutMatrix_01(){
	Matrix E(3,3);
	E(1,1)=0.999999998306989;E(1,2)=-5.33807310308936e-05;E(1,3)=-2.31628936941766e-05;
	E(2,1)=5.33814884667738e-05;E(2,2)=0.99999999804053;E(2,3)=3.27010136421668e-05;
	E(3,1)=2.31611480447759e-05;E(3,2)=-3.27022500565821e-05;E(3,3)=0.999999999197062;

	_assert(m_equals(E,NutMatrix(2025), 1e-10));

	return 0;
}

int PoleMatrix_01(){
	Matrix E(3,3);
	E(1,1)=0.54030230586814;E(1,2)=0;E(1,3)=0.841470984807897;
	E(2,1)=0;E(2,2)=1;E(2,3)=0;
	E(3,1)=-0.841470984807897;E(3,2)=0;E(3,3)=0.54030230586814;

	_assert(m_equals(E,PoleMatrix(1,0), 1e-10));

	return 0;
}

int PrecMatrix_01(){
	Matrix E(3,3);
	E(1,1)=0.999999999999778;E(1,2)=-6.11727857776369e-07;E(1,3)=-2.66195216523163e-07;
	E(2,1)=6.11727857776369e-07;E(2,2)=0.999999999999813;E(2,3)=-8.14195151602225e-14;
	E(3,1)=2.66195216523163e-07;E(3,2)=-8.1419514393827e-14;E(3,3)=0.999999999999965;

	_assert(m_equals(E,PrecMatrix(2024,2025), 1e-10));

	return 0;
}

int gmst_01(){
	_assert(fabs(gmst(2025)-4.39293397920984)<1e-10);

	return 0;
}

int gast_01(){
	_assert(fabs(gast(2025)-4.392987359940904)<1e-10);

	return 0;
}

int MeasUpdate_01(){
	Matrix x(3,1);
	x(1,1)=1;x(2,1)=1;x(3,1)=1;

	Matrix G(3);
	G(1)=2;G(2)=2;G(2)=2;

	Matrix P(3,3);
	P(1,1)=1;P(1,2)=2;P(1,3)=3;
	P(2,1)=4;P(2,2)=5;P(2,3)=6;
	P(3,1)=7;P(3,2)=8;P(3,3)=9;

	auto[A,B,C]=MeasUpdate(x,1.5,0.5,0.005,G,P,3);

	Matrix A1(3);
	A1(1)=0.066666657407409;A1(2)=0.166666643518522;A1(3)=0.266666629629635;
	Matrix B1(3);
	B1(1)=1.066666657407409;B1(2)=1.166666643518522;B1(3)=1.266666629629635;
	Matrix C1(3,3);
	C1(1,1)=-0.599999777777809;C1(1,2)=0.000000277777739;C1(1,3)=0.600000333333287;
	C1(2,1)=0.000000555555479;C1(2,2)=0.000000694444348;C1(2,3)=0.000000833333218;
	C1(3,1)=0.600000888888765;C1(3,2)=0.000001111110957;C1(3,3)=-0.599998666666853;

	_assert(m_equals(A,transpose(A1), 1e-10));
	_assert(m_equals(B,transpose(B1), 1e-10));
	_assert(m_equals(C,C1, 1e-10));

	return 0;
}

int all_tests()
{
    _verify(m_sum_01);
    _verify(m_sub_01);
    _verify(m_zeros_01);
    _verify(m_zeros_02);
	_verify(m_mul_01);
	_verify(m_div_01);
	_verify(m_eq_01);
	_verify(m_eye_01);
	_verify(m_trans_01);
	_verify(m_inv_01);
	_verify(m_sum_02);
	_verify(m_sub_02);
	_verify(m_mul_02);
	_verify(m_div_02);
	_verify(m_norm_01);
	_verify(m_dot_01);
	_verify(m_cross_01);
	_verify(m_extractVector_01);
	_verify(m_extractRow_01);
	_verify(m_extractColumn_01);
	_verify(m_unionVector_01);
	_verify(m_assignRow_01);
	_verify(m_assignColumn_01);
	_verify(R_x_01);
	_verify(R_y_01);
	_verify(R_z_01);
	_verify(AccelPointMass_01);
	_verify(Cheb3D_01);
	_verify(EccAnom_01);
	_verify(Frac_01);
	_verify(MeanObliquity_01);
	_verify(Mjday_01);
	_verify(Mjday_02);
	_verify(MjdayTDB_01);
	_verify(Position_01);
	_verify(sign_01);
	_verify(AzElPa_01);
	_verify(IERS_01);
	_verify(IERS_02);
	_verify(Legendre_01);
	_verify(NutAngles_01);
	_verify(TimeUpdate_01);
	_verify(TimeUpdate_02);
	_verify(JPL_Eph_DE430_01);
	_verify(AccelHarmonic_01);
	_verify(EqnEquinox_01);
	_verify(LTC_01);
	_verify(NutMatrix_01);
	_verify(PoleMatrix_01);
	_verify(PrecMatrix_01);
	_verify(gmst_01);
	_verify(gast_01);
	_verify(MeanObliquity_01);

    return 0;
}


int main()
{
	eop19620101(21413);
	DE430Coeff(2285,1020);
	GGM03S(181);
	
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
