#include "..\include\matrix.h"

Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}

double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}

Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}

ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix& Matrix::operator * (Matrix &m) {
	if (this->n_column != m.n_row) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, m.n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
		for(int k = 1; k <= m.n_column; k++){
			(*m_aux)(i,k) = 0;
			for(int j = 1; j <= this->n_column; j++) {
				(*m_aux)(i,k) += (*this)(i,j) * m(j,k);
			}
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator = (Matrix &m) {
	if (this == &m) return *this;

	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
		
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*this)(i,j) = m(i,j);
		}
	}
	
	return *this;
}

Matrix& Matrix::operator / (Matrix &m) {
	return (*this) * inv(m) ;
}

Matrix& eye(const int n_row) {
	Matrix *m_aux = new Matrix(n_row, n_row);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_row; j++) {
			if(i==j){(*m_aux)(i,j) = 1;}
			else {(*m_aux)(i,j) = 0;}
		}
	}

	return (*m_aux);
}

Matrix& transpose(Matrix &m) {
	Matrix *m_aux = new Matrix(m.n_column,m.n_row);
	
	for(int i = 1; i <= m.n_row; i++) {
		for(int j = 1; j <= m.n_column; j++){
			(*m_aux)(j,i) = m(i,j);
		}
	}
	
	return (*m_aux);
}

Matrix& inv(Matrix &m) {
	if (m.n_row != m.n_column) {
        cout << "Matrix inversion error: Matrix is not square\n";
        exit(EXIT_FAILURE);
    }

	int n=m.n_row;

	Matrix *m_aux = new Matrix(n,2*n);
	
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++){
			(*m_aux)(i,j) = m(i,j);
		}
		for(int j = n+1; j <= 2*n; j++){
			(*m_aux)(i,j) = (i==(j-n)) ? 1.0 : 0.0;
		}
	}

	for (int i = 1; i <= n; i++) {
        // Partial Pivoting: Find the row with the maximum element in this column
        int pivotRow = i;
        for (int k = i + 1; k <= n; k++) {
            if (fabs((*m_aux)(k, i)) > fabs((*m_aux)(pivotRow, i))) {
                pivotRow = k;
            }
        }

        // Swap rows if necessary
        if (pivotRow != i) {
            for (int j = 1; j <= 2 * n; j++) {
                swap((*m_aux)(i, j), (*m_aux)(pivotRow, j));
            }
        }

        // Check for singular matrix (non-invertible)
        if (fabs((*m_aux)(i, i)) < 1e-10) {
            cout << "Matrix inversion error: Matrix is singular (non-invertible)\n";
            exit(EXIT_FAILURE);
        }

        // Normalize the pivot row (make leading coefficient = 1)
        double pivotValue = (*m_aux)(i, i);
        for (int j = 1; j <= 2 * n; j++) {
            (*m_aux)(i, j) /= pivotValue;
        }

        // Eliminate all other elements in the current column
        for (int k = 1; k <= n; k++) {
            if (k != i) {
                double factor = (*m_aux)(k, i);
                for (int j = 1; j <= 2 * n; j++) {
                    (*m_aux)(k, j) -= factor * (*m_aux)(i, j);
                }
            }
        }
    }

	Matrix *inverse = new Matrix(n, n);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*inverse)(i, j) = (*m_aux)(i, j + n);
        }
    }
	
	return (*inverse);
}