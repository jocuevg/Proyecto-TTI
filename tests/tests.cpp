#include "..\include\matrix.h"
#include <cstdio>
#include <cmath>

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

int all_tests()
{
    _verify(m_sum_01);
    _verify(m_sub_01);
    _verify(m_zeros_01);
    _verify(m_zeros_02);
	_verify(m_mul_01);
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

    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
