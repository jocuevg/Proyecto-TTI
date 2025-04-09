//$Source$
//------------------------------------------------------------------------------
// matrix
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file matrix.cpp
* @brief Programacion de operaciones de matrices.
*
* @author Jose Cuevas Gil de Gomez
* @bug No hay.
*/
//------------------------------------------------------------------------------
#include "..\include\matrix.hpp"

Matrix::Matrix(const int n_row, const int n_column)
{
	if (n_row <= 0 || n_column <= 0)
	{
		cout << "Matrix create: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **)malloc(n_row * sizeof(double *));

	if (this->data == NULL)
	{
		cout << "Matrix create: error in data\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < n_row; i++)
	{
		this->data[i] = (double *)malloc(n_column * sizeof(double));
	}
}

Matrix::Matrix(const int v_size)
{
	if (v_size < 0)
	{
		cout << "Vector create: error in n_size\n";
		exit(EXIT_FAILURE);
	}

	this->n_row = 1;
	this->n_column = v_size;
	this->data = (double **)malloc(n_row * sizeof(double *));

	if (this->data == NULL)
	{
		cout << "Vector create: error in data\n";
		exit(EXIT_FAILURE);
	}

	this->data[0] = (double *)calloc(v_size, sizeof(double));
}

double &Matrix::operator()(const int row, const int column)
{
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column)
	{
		cout << "Matrix get: error in row/column\n";
		exit(EXIT_FAILURE);
	}

	return this->data[row - 1][column - 1];
}

double &Matrix::operator()(const int n)
{
	if (n <= 0 || n > this->n_row * this->n_column)
	{
		cout << "Vector get: error in row/column\n";
		exit(EXIT_FAILURE);
	}

	return this->data[(n - 1) / this->n_column][(n - 1) % this->n_column];
}

Matrix &Matrix::operator+(Matrix &m)
{
	if (this->n_row != m.n_row || this->n_column != m.n_column)
	{
		cout << "Matrix sum: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) + m(i, j);
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator-(Matrix &m)
{
	if (this->n_row != m.n_row || this->n_column != m.n_column)
	{
		cout << "Matrix sub: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) - m(i, j);
		}
	}

	return *m_aux;
}

ostream &operator<<(ostream &o, Matrix &m)
{
	for (int i = 1; i <= m.n_row; i++)
	{
		for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i, j));
		o << "\n";
	}

	return o;
}

Matrix &zeros(const int n_row, const int n_column)
{
	Matrix *m_aux = new Matrix(n_row, n_column);

	for (int i = 1; i <= n_row; i++)
	{
		for (int j = 1; j <= n_column; j++)
		{
			(*m_aux)(i, j) = 0;
		}
	}

	return (*m_aux);
}

Matrix &zeros(const int n)
{
	Matrix *m_aux = new Matrix(n);

	for (int i = 1; i <= n; i++)
	{
		(*m_aux)(i) = 0;
	}

	return (*m_aux);
}

Matrix &Matrix::operator*(Matrix &m)
{
	if (this->n_column != m.n_row)
	{
		cout << "Matrix sub: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(this->n_row, m.n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int k = 1; k <= m.n_column; k++)
		{
			(*m_aux)(i, k) = 0;
			for (int j = 1; j <= this->n_column; j++)
			{
				(*m_aux)(i, k) += (*this)(i, j) * m(j, k);
			}
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator=(Matrix &m)
{
	if (this == &m)
		return *this;

	if (this->n_row != m.n_row || this->n_column != m.n_column)
	{
		cout << "Matrix sub: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*this)(i, j) = m(i, j);
		}
	}

	return *this;
}

Matrix &Matrix::operator/(Matrix &m)
{
	return (*this) * inv(m);
}

Matrix &eye(const int n_row)
{
	Matrix *m_aux = new Matrix(n_row, n_row);

	for (int i = 1; i <= n_row; i++)
	{
		for (int j = 1; j <= n_row; j++)
		{
			if (i == j)
			{
				(*m_aux)(i, j) = 1;
			}
			else
			{
				(*m_aux)(i, j) = 0;
			}
		}
	}

	return (*m_aux);
}

Matrix &transpose(Matrix &m)
{
	Matrix *m_aux = new Matrix(m.n_column, m.n_row);

	for (int i = 1; i <= m.n_row; i++)
	{
		for (int j = 1; j <= m.n_column; j++)
		{
			(*m_aux)(j, i) = m(i, j);
		}
	}

	return (*m_aux);
}

Matrix &inv(Matrix &m)
{
	if (m.n_row != m.n_column)
	{
		cout << "Matrix inversion error: Matrix is not square\n";
		exit(EXIT_FAILURE);
	}

	int n = m.n_row;

	Matrix *m_aux = new Matrix(n, 2 * n);

	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			(*m_aux)(i, j) = m(i, j);
		}
		for (int j = n + 1; j <= 2 * n; j++)
		{
			(*m_aux)(i, j) = (i == (j - n)) ? 1.0 : 0.0;
		}
	}

	for (int i = 1; i <= n; i++)
	{
		// Partial Pivoting: Find the row with the maximum element in this column
		int pivotRow = i;
		for (int k = i + 1; k <= n; k++)
		{
			if (fabs((*m_aux)(k, i)) > fabs((*m_aux)(pivotRow, i)))
			{
				pivotRow = k;
			}
		}

		// Swap rows if necessary
		if (pivotRow != i)
		{
			for (int j = 1; j <= 2 * n; j++)
			{
				swap((*m_aux)(i, j), (*m_aux)(pivotRow, j));
			}
		}

		// Check for singular matrix (non-invertible)
		if (fabs((*m_aux)(i, i)) < 1e-10)
		{
			cout << "Matrix inversion error: Matrix is singular (non-invertible)\n";
			exit(EXIT_FAILURE);
		}

		// Normalize the pivot row (make leading coefficient = 1)
		double pivotValue = (*m_aux)(i, i);
		for (int j = 1; j <= 2 * n; j++)
		{
			(*m_aux)(i, j) /= pivotValue;
		}

		// Eliminate all other elements in the current column
		for (int k = 1; k <= n; k++)
		{
			if (k != i)
			{
				double factor = (*m_aux)(k, i);
				for (int j = 1; j <= 2 * n; j++)
				{
					(*m_aux)(k, j) -= factor * (*m_aux)(i, j);
				}
			}
		}
	}

	Matrix *inverse = new Matrix(n, n);
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			(*inverse)(i, j) = (*m_aux)(i, j + n);
		}
	}

	return (*inverse);
}

Matrix &Matrix::operator+(const double x)
{
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) + x;
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator-(const double x)
{
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) - x;
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator*(const double x)
{
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) * x;
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator/(const double x)
{
	if (x == 0)
	{
		cout << "Matrix div: error in x is zero\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) / x;
		}
	}

	return *m_aux;
}

double norm(Matrix &m)
{

	if (m.n_row != 1)
	{
		cout << "Matrix inversion error: Matrix is not square\n";
		exit(EXIT_FAILURE);
	}
	double res = 0;

	for (int i = 1; i <= m.n_column; i++)
	{
		res += pow(m(1, i), 2);
	}

	return sqrt(res);
}

double dot(Matrix &m, Matrix &n)
{

	if (m.n_row != 1 || n.n_row != 1 || m.n_column != n.n_column)
	{
		cout << "Matrix inversion error: Matrix is not square\n";
		exit(EXIT_FAILURE);
	}
	double res = 0;

	for (int i = 1; i <= m.n_column; i++)
	{
		res += m(1, i) * n(1, i);
	}

	return res;
}

Matrix &cross(Matrix &m, Matrix &n)
{

	if (m.n_row != 1 || n.n_row != 1 || m.n_column != n.n_column || m.n_column != 3)
	{
		cout << "Matrix inversion error: Matrix is not square\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(3);

	(*m_aux)(1, 1) = m(1, 2) * n(1, 3) - m(1, 3) * n(1, 2);
	(*m_aux)(1, 2) = m(1, 3) * n(1, 1) - m(1, 1) * n(1, 3);
	(*m_aux)(1, 3) = m(1, 1) * n(1, 2) - m(1, 2) * n(1, 1);

	return *m_aux;
}

Matrix &extract_vector(Matrix &m, const int n, const int k)
{

	if (m.n_row != 1 || n <= 0 || k < n || k > m.n_column)
	{
		cout << "Vector extract: error in indexes\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(k - n + 1);

	for (int i = n; i <= k; i++)
	{
		(*m_aux)(i - n + 1) = m(i);
	}

	return *m_aux;
}

Matrix &extract_row(Matrix &m, const int n)
{
	if (n <= 0 || n > m.n_row)
	{
		cout << "Matrix get: error in row/column\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(1, m.n_column);

	for (int j = 1; j <= m.n_column; j++)
	{
		(*m_aux)(1, j) = m(n, j);
	}

	return (*m_aux);
}

Matrix &extract_column(Matrix &m, const int n)
{
	if (n <= 0 || n > m.n_column)
	{
		cout << "Matrix get: error in row/column\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(1, m.n_row);

	for (int j = 1; j <= m.n_row; j++)
	{
		(*m_aux)(1, j) = m(j, n);
	}

	return (*m_aux);
}

Matrix &union_vector(Matrix &m, Matrix &k)
{
	if (m.n_row != 1 || k.n_row != 1)
	{
		cout << "Vector union: error in vector\n";
		exit(EXIT_FAILURE);
	}

	int t=m.n_column + k.n_column;

	Matrix *m_aux = new Matrix(t);

	for (int i = 1; i <= m.n_column; i++)
	{
		(*m_aux)(i) = m(i);
	}
	for (int i = 1; i <= k.n_column; i++)
	{
		(*m_aux)(i+m.n_column) = k(i);
	}

	return *m_aux;
}

Matrix & assign_column(Matrix&m,Matrix&k, const int n){
	if (k.n_row != 1 || k.n_column != m.n_row || n<=0 || n>m.n_column)
	{
		cout << "Assign column: error in vector/index\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 1; i <= m.n_row; i++)
	{
		m(i,n)=k(i);
	}

	return m;
}

Matrix & assign_row(Matrix&m,Matrix&k,const int n){
	if (k.n_row != 1 || k.n_column != m.n_column || n<=0 || n>m.n_row)
	{
		cout << "Assign column: error in vector/index\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 1; i <= m.n_column; i++)
	{
		m(n,i)=k(i);
	}

	return m;
}