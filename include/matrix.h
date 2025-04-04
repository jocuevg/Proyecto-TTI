#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix {
	private:
	double **data;


public:
    int n_row, n_column;

    // Parameterized constructor
    Matrix(const int v_size);
    Matrix(const int n_row, const int n_column);
	
	// Member operators
	double& operator () (const int n);
	double& operator () (const int row, const int column);
	Matrix& operator + (Matrix &m);
	Matrix& operator - (Matrix &m);
	Matrix& operator * (Matrix &m);
	Matrix& operator / (Matrix &m);
	Matrix& operator = (Matrix &m);
	Matrix& operator + (const double x);
	Matrix& operator - (const double x);
	Matrix& operator * (const double x);
	Matrix& operator / (const double x);
	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int n_row, const int n_column);
Matrix& zeros(const int n);
Matrix& eye(const int n_row);
Matrix& transpose(Matrix &m);
Matrix& inv(Matrix &m);
double norm(Matrix &m);
double dot(Matrix &m,Matrix &n);
Matrix & cross(Matrix &m,Matrix &n);
Matrix & extract_vector(Matrix&m,const int n, const int k);
Matrix & extract_row(Matrix&m,const int n);
Matrix & extract_column(Matrix&m,const int n);
Matrix & union_vector(Matrix&m,Matrix&k,const int n);
Matrix & assign_column(Matrix&m,const int n);
Matrix & assign_row(Matrix&m,const int n);
#endif