/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#include "Matrix.h"

Matrix::Matrix()
{
	data = NULL;
}

Matrix::Matrix(const int m, const int n)
{
	data = gsl_matrix_calloc(m,n);
}

Matrix::Matrix(const int n)
{
	data = gsl_matrix_calloc(n,n);
}

Matrix::Matrix(const Matrix & A)
{
	if(A.data != NULL)
	{
		data = gsl_matrix_alloc(A.data->size1, A.data->size2);
		gsl_matrix_memcpy(data, A.data);
	}
	else
	{
		data = NULL;
	}
}

Matrix::~Matrix()
{
	if(data != NULL)
	{
		gsl_matrix_free(data);
	}
}

double Matrix::get(const int i, const int j) const
{
	return gsl_matrix_get(data, i, j);
}

void Matrix::set(const double v, const int i, const int j)
{
	gsl_matrix_set(data, i, j, v);
}

int Matrix::rows() const
{
	return data->size1;
}

int Matrix::cols() const
{
	return data->size2;
}

void Matrix::row(RowVector & result, const int i) const
{
	for(int j=0; j<data->size2; ++j)
	{
		result.set(gsl_matrix_get(data, i, j), j);
	}
}

void Matrix::sortColumns()
{
	int m = data->size1;
	int n = data->size2;

	double *sizes = new double[n];

	// Compute the sizes of the columns
	for(int j=0; j<n; ++j)
	{
		double size = 0;
		for(int i=0; i<m; ++i)
		{
			double tmp = gsl_matrix_get(data, i, j);
			tmp *= tmp;
			size += tmp;
		}
		sizes[j] = size;
	}

	// Selection sort
	double tmp;
	int iMax;

	for(int i=0; i<n-1; ++i)
	{
		iMax = i;
		for(int j=i+1; j<n; ++j)
		{
			if(sizes[j] > sizes[iMax])
			{
				iMax = j;
			}
		}

		//Exchange the columns
		if(iMax != i)
		{
			gsl_matrix_swap_columns(data, i, iMax);
			tmp = sizes[i];
			sizes[i] = sizes[iMax];
			sizes[iMax] = tmp;
		}
	}

	delete[] sizes;
}

int Matrix::rank() const
{
	int m = data->size1;
	int n = data->size2;

	Matrix temp(*this);

	gsl_vector *work = gsl_vector_alloc(n);
	gsl_vector *S = gsl_vector_alloc(n);
	gsl_matrix *X = gsl_matrix_alloc(n,n);
	gsl_matrix *V = gsl_matrix_alloc(n,n);

	gsl_linalg_SV_decomp_mod(temp.data, X, V, S, work);

	int r = 0;
	double tmp;
	for(int i=0; i<n; ++i)
	{
		tmp = gsl_vector_get(S, i);
		if(tmp < THRESHOLD_LOW)
			break;
		else
			++r;
	}

	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(X);
	gsl_matrix_free(V);

	return r;
}

void Matrix::transpose(Matrix & result) const
{
	gsl_matrix_transpose_memcpy(result.data, data);
}

void Matrix::svd(Matrix & U) const
{
	U = *this;
	int m = data->size1;
	int n = data->size2;

	gsl_vector *work = gsl_vector_alloc(n);
	gsl_vector *S = gsl_vector_alloc(n);
	gsl_matrix *V = gsl_matrix_alloc(n,n);

	gsl_linalg_SV_decomp(U.data, V, S, work);

	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(V);
}

void Matrix::neg(Matrix & result) const
{
	gsl_matrix_memcpy(result.data, data);
	gsl_matrix_scale(result.data, -1.0);
}

void Matrix::neg_assign()
{
	gsl_matrix_scale(data, -1.0);
}

void Matrix::inverse(Matrix & result) const
{
	// only square matrix
	int m = data->size1;
	int n = data->size2;

	if(m != n)
	{
		printf("Not a square matrix.\n");
		return;
	}

	// We use GSL library.
	gsl_matrix *A = gsl_matrix_alloc(m,m);
	gsl_permutation *p = gsl_permutation_alloc(m);
	gsl_matrix *invA = gsl_matrix_alloc(m,m);

	// Make a copy the matrix.
	gsl_matrix_memcpy(A, data);

	int *signum = new int[m];

	gsl_linalg_LU_decomp(A, p, signum);
	gsl_linalg_LU_invert(A, p, invA);

	gsl_matrix_memcpy(result.data, invA);

	gsl_matrix_free(A);
	gsl_permutation_free(p);
	gsl_matrix_free(invA);
	delete signum;
}

void Matrix::inverse_assign()
{
	// only square matrix
	int m = data->size1;
	int n = data->size2;

	if(m != n)
	{
		printf("Not a square matrix.\n");
		return;
	}

	// We use GSL library.
	gsl_matrix *A = gsl_matrix_alloc(m,m);
	gsl_permutation *p = gsl_permutation_alloc(m);
	gsl_matrix *invA = gsl_matrix_alloc(m,m);

	// Make a copy the matrix.
	gsl_matrix_memcpy(A, data);

	int *signum = new int[m];

	gsl_linalg_LU_decomp(A, p, signum);
	gsl_linalg_LU_invert(A, p, invA);

	gsl_matrix_memcpy(data, invA);

	gsl_matrix_free(A);
	gsl_permutation_free(p);
	gsl_matrix_free(invA);
	delete signum;
}

Matrix & Matrix::operator += (const Matrix & A)
{
	gsl_matrix_add(data, A.data);

	return *this;
}

Matrix & Matrix::operator -= (const Matrix & A)
{
	gsl_matrix_sub(data, A.data);

	return *this;
}

Matrix & Matrix::operator *= (const Matrix & A)
{
	int m = data->size1;
	int n = A.data->size2;
	int k = data->size2;

	Matrix result(m,n);

	for(int i=0; i<m; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			double tmp = 0;
			for(int p=0; p<k; ++p)
			{
				tmp += gsl_matrix_get(data, i, p) * gsl_matrix_get(A.data, p, j);
			}
			gsl_matrix_set(result.data, i, j, tmp);
		}
	}

	*this = result;
	return *this;
}

Matrix Matrix::operator + (const Matrix & A) const
{
	Matrix result = *this;

	result += A;

	return result;
}

Matrix Matrix::operator - (const Matrix & A) const
{
	Matrix result = *this;

	result -= A;

	return result;
}

Matrix Matrix::operator * (const Matrix & A) const
{
	int m = data->size1;
	int n = A.data->size2;
	int k = data->size2;

	Matrix result(m,n);

	for(int i=0; i<m; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			double tmp = 0;
			for(int p=0; p<k; ++p)
			{
				tmp += gsl_matrix_get(data, i, p) * gsl_matrix_get(A.data, p, j);
			}
			gsl_matrix_set(result.data, i, j, tmp);
		}
	}

	return result;
}

void Matrix::QR(Matrix & D)
{
	int m = data->size1;
	int n = data->size2;

	int i,j,k,l;

	k = 0;
	double s,t,x,r;

	for(l=0; l<n; ++l)
	{
		if(k == m)
		{
			gsl_matrix_set(D.data, 0, l, gsl_matrix_get(data, k, l));
			break;
		}
		s = 0;

		for (i=k; i<m; ++i)
		{
			x = gsl_matrix_get(data, i, l);
			s += x*x;
		}
		s = sqrt(s);

		if(s == 0)
		{
			gsl_matrix_set(D.data, 0, 0, 0);
			continue;
		}

		t = gsl_matrix_get(data, k, l);
		r = 1/sqrt( s*(s+fabs(t)) );
		if(t < 0)
			s = -s;

		gsl_matrix_set(D.data, 0, l, -s);
		gsl_matrix_set(data, k, k, r * (t + s));

		for(i=k+1; i<m; ++i)
		{
			double tmp = gsl_matrix_get(data, i, l) * r;
			gsl_matrix_set(data, i, k, tmp);
		}

		for(j=l+1; j<n; ++j)
		{
			t = 0;
			for(i=k; i<m; ++i)
				t += gsl_matrix_get(data, i, k) * gsl_matrix_get(data, i, j);

			for(i=k; i<m; ++i)
			{
				double tmp = gsl_matrix_get(data, i, j);
				gsl_matrix_set(data, i, j, tmp - gsl_matrix_get(data, i, k)*t);
			}
		}

		++k;
    }
}

void Matrix::QRfactor(Matrix & Q)
{
	int m = data->size1;
	int n = data->size2;

	if(n == 1)
    {
		gsl_matrix_set(Q.data, 0, 0, 1);
		return;
    }

	Q = *this;

	Matrix D(1,n);

	Q.QR(D);

	Matrix V(1,n);

	gsl_matrix_set(V.data, 0, n-2, gsl_matrix_get(Q.data, n-2, n-2));
	gsl_matrix_set(V.data, 0, n-1, gsl_matrix_get(Q.data, n-1, n-2));

	gsl_matrix_set(Q.data, n-1, n-2, - gsl_matrix_get(V.data, 0, n-2) * gsl_matrix_get(V.data, 0, n-1));
	gsl_matrix_set(Q.data, n-2, n-1, gsl_matrix_get(Q.data, n-1, n-2));

	gsl_matrix_set(Q.data, n-1, n-1, 1.0 - gsl_matrix_get(V.data, 0, n-1) * gsl_matrix_get(V.data, 0, n-1));
	gsl_matrix_set(Q.data, n-2, n-2, 1.0 - gsl_matrix_get(V.data, 0, n-2) * gsl_matrix_get(V.data, 0, n-2));

	int k, row, col, i;
	double a, b;
	for(k = n-3; k>=0; --k)
    {
		for(row = k; row<n; ++row)
			gsl_matrix_set(V.data, 0, row, gsl_matrix_get(Q.data, row, k));

		for(col = k+1; col<n; ++col)
		{
			a = 0;
			for(i = k+1; i<n; ++i)
				a += gsl_matrix_get(V.data, 0, i) * gsl_matrix_get(Q.data, i, col);

			for(row = k+1; row<n; ++row)
			{
				b = gsl_matrix_get(Q.data, row, col);
				gsl_matrix_set(Q.data, row, col, b - a * gsl_matrix_get(V.data, 0, row));
			}

			gsl_matrix_set(Q.data, k, col, -a * gsl_matrix_get(V.data, 0, k));
		}

		for(i = k+1; i<n; ++i)
			gsl_matrix_set(Q.data, i, k, - gsl_matrix_get(V.data, 0, k) * gsl_matrix_get(V.data, 0, i));

		gsl_matrix_set(Q.data, k, k, 1.0 - gsl_matrix_get(V.data, 0, k) * gsl_matrix_get(V.data, 0, k));
    }
}

void Matrix::output(FILE *fp) const
{
	int m = data->size1;
	int n = data->size2;

	fprintf(fp, "==========\n");
	for(int i=0; i<m; ++i)
	{
		for(int j=0; j<n-1; ++j)
		{
			fprintf(fp, "%lf, ", gsl_matrix_get(data, i, j));
		}
		fprintf(fp, "%lf\n", gsl_matrix_get(data, i, n-1));
	}
	fprintf(fp, "==========\n");
}

Matrix & Matrix::operator = (const Matrix & A)
{
	if(A.data == NULL)
	{
		if(data != NULL)
		{
			gsl_matrix_free(data);
			data = NULL;
		}
	}
	else
	{
		if(data != NULL)
			gsl_matrix_free(data);

		data = gsl_matrix_alloc(A.data->size1, A.data->size2);

		gsl_matrix_memcpy(data, A.data);
	}
	return *this;
}










// class RowVector

RowVector::RowVector()
{
}

RowVector::RowVector(const int n)
{
	Matrix mat(1,n);
	vec = mat;
}

RowVector::RowVector(const RowVector & v)
{
	vec = v.vec;
}

RowVector::~RowVector()
{
}

void RowVector::set(const double v, const int pos)
{
	vec.set(v, 0, pos);
}

double RowVector::get(const int pos) const
{
	return vec.get(0, pos);
}

int RowVector::size() const
{
	return vec.cols();
}

void RowVector::transpose(ColVector & result) const
{
	Matrix mat(vec.cols(), 1);
	vec.transpose(mat);
	result.vec = mat;
}

void RowVector::neg(RowVector & result) const
{
	result = *this;
	result.vec.neg_assign();
}

void RowVector::neg_assign()
{
	vec.neg_assign();
}

void RowVector::dump(FILE *fp) const
{
	fprintf(fp, "[ ");
	for(int i=0; i<vec.cols()-1; ++i)
	{
		fprintf(fp, "%lf, ", get(i));
	}
	fprintf(fp, "%lf ]\n", get(vec.cols()-1));
}

double RowVector::innerProd(const RowVector & v) const
{
	int n = vec.cols();

	if(n != v.size())
	{
		printf("Vector dimensions do not coincide.\n");
		return INVALID;
	}

	double result = 0;
	for(int i=0; i<n; ++i)
	{
		result += get(i)*v.get(i);
	}

	return result;
}

double RowVector::EuclideanNorm() const
{
	int n = vec.cols();
	double result = 0;

	for(int i=0; i<n; ++i)
	{
		result += get(i)*get(i);
	}

	return sqrt(result);
}

void RowVector::normalize()
{
	double norm = EuclideanNorm();
	int n = vec.cols();

	for(int i=0; i<n; ++i)
	{
		double tmp = get(i) / norm;
		set(tmp, i);
	}
}

bool RowVector::operator == (const RowVector & v) const
{
	if(vec.cols() == v.size())
	{
		for(int i=0; i<vec.cols(); ++i)
		{
			double d = vec.get(0,i) - v.get(i);
			if(!(d <= THRESHOLD_LOW && d >= -THRESHOLD_LOW))
			{
				return false;
			}
		}

		return true;
	}
	else
	{
		return false;
	}
}

RowVector & RowVector::operator += (const RowVector & v)
{
	vec += v.vec;
	return *this;
}

RowVector & RowVector::operator -= (const RowVector & v)
{
	vec -= v.vec;
	return *this;
}

RowVector RowVector::operator + (const RowVector & v) const
{
	RowVector result = *this;
	result += v;
	return result;
}

RowVector RowVector::operator - (const RowVector & v) const
{
	RowVector result = *this;
	result -= v;
	return result;
}

RowVector & RowVector::operator = (const RowVector & v)
{
	if(this == &v)
		return *this;

	vec = v.vec;
	return *this;
}









// class ColVector

ColVector::ColVector()
{
}

ColVector::ColVector(const int n)
{
	Matrix mat(n,1);
	vec = mat;
}

ColVector::ColVector(const ColVector & v)
{
	vec = v.vec;
}

ColVector::~ColVector()
{
}

void ColVector::set(const double v, const int pos)
{
	vec.set(v, pos, 0);
}

double ColVector::get(const int pos) const
{
	return vec.get(pos, 0);
}

int ColVector::size() const
{
	return vec.rows();
}

void ColVector::transpose(RowVector & result) const
{
	Matrix mat(1, vec.cols());
	vec.transpose(mat);
	result.vec = mat;
}

void ColVector::neg(ColVector & result) const
{
	result = *this;
	result.vec.neg_assign();
}

void ColVector::neg_assign()
{
	vec.neg_assign();
}

void ColVector::mul(ColVector & result, const Matrix & m) const
{
	int rows = m.rows();
	int cols = m.cols();

	if(cols != vec.rows())
	{
		printf("Vector multiplication error: invalid dimensions.\n");
		return;
	}

	for(int i=0; i<rows; ++i)
	{
		double sum = 0;
		for(int j=0; j<cols; ++j)
		{
			sum += m.get(i,j) * vec.get(j,0);
		}

		result.vec.set(sum, i, 0);
	}
}

void ColVector::mul_assign(const Matrix & m)
{
	ColVector result(m.rows());

	mul(result, m);
	*this = result;
}

ColVector & ColVector::operator += (const ColVector & v)
{
	vec += v.vec;
	return *this;
}

ColVector & ColVector::operator -= (const ColVector & v)
{
	vec -= v.vec;
	return *this;
}

ColVector ColVector::operator + (const ColVector & v) const
{
	ColVector result = *this;
	result += v;
	return result;
}

ColVector ColVector::operator - (const ColVector & v) const
{
	ColVector result = *this;
	result -= v;
	return result;
}

ColVector & ColVector::operator = (const ColVector & v)
{
	if(this == &v)
		return *this;

	vec = v.vec;
	return *this;
}
