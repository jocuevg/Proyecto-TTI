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
	
	Matrix A(f);
	A(1,1) = 1; A(1,2) =  2; A(1,3) = 3; 
	
	Matrix B(f);
	B(1,1) = 4; B(1,2) =  5; B(1,3) = 6;
	
	Matrix C(f);
	C(1,1) = 0.0463899400720928; C(1,2) =  0.0419499176126264; C(1,3) = 0.0375098951531599; 
	
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

	Matrix A(3);
	A(1)=1;A(2)=1;A(3)=1;
	tuple<double,double,Matrix,Matrix> R = AzElPa(A);
	
	Matrix B(3);
	B(1)=0.5;B(2)=-0.5;B(3)=0;
	Matrix C(3);
	C(1)=-0.23570226039;C(2)=-0.23570226039;C(3)=0.47140452079;

	tuple<double,double,Matrix,Matrix> D = make_tuple(0.78539816339,0.61547970867,B,C);
	
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

    return 0;
}


int main()
{
	eop19620101(21413);
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
