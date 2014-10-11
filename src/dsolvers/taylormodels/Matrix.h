//modified by Sicun Gao
/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#ifndef MATRIX_H_
#define MATRIX_H_

#include "flowstar_include.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

// The matrix class is implemented based on the data structure of gsl matrix.

class RowVector;
class ColVector;

class Matrix
{
private:
	gsl_matrix *data;
public:
	Matrix();
	Matrix(const int m, const int n);	// Create an m x n matrix, all of the entries are 0.
	Matrix(const int n);				// Create an n x n matrix, all of the entries are 0.
	Matrix(const Matrix & A);
	~Matrix();

	double get(const int i, const int j) const;				// Get the entry at position [i,j].
	void set(const double v, const int i, const int j);		// Set A[i,j] = v.

	int rows() const;
	int cols() const;

	void row(RowVector & result, const int i) const;			// Return the (i+1)-st row.

	void sortColumns();		// Sort the columns by size in descending order.
	int rank() const;

	void neg(Matrix & result) const;
	void neg_assign();

	void inverse(Matrix & result) const;
	void inverse_assign();

	void transpose(Matrix & result) const;
	void svd(Matrix & U) const;

	void QR(Matrix & D);
	void QRfactor(Matrix & Q);

	void output(FILE *fp) const;

	Matrix & operator += (const Matrix & A);
	Matrix & operator -= (const Matrix & A);
	Matrix & operator *= (const Matrix & A);

	Matrix operator + (const Matrix & A) const;
	Matrix operator - (const Matrix & A) const;
	Matrix operator * (const Matrix & A) const;

	Matrix & operator = (const Matrix & A);
};

class RowVector
{
private:
	Matrix vec;
public:
	RowVector();
	RowVector(const int n);
	RowVector(const RowVector & v);
	~RowVector();

	void set(const double v, const int pos);
	double get(const int pos) const;
	int size() const;

	void transpose(ColVector & result) const;

	void neg(RowVector & result) const;
	void neg_assign();

	void dump(FILE *fp) const;

	double innerProd(const RowVector & v) const;
	double EuclideanNorm() const;
	void normalize();

	bool operator == (const RowVector & v) const;

	RowVector & operator += (const RowVector & v);
	RowVector & operator -= (const RowVector & v);
	RowVector operator + (const RowVector & v) const;
	RowVector operator - (const RowVector & v) const;

	RowVector & operator = (const RowVector & v);

	friend class ColVector;
};

class ColVector
{
private:
	Matrix vec;
public:
	ColVector();
	ColVector(const int n);
	ColVector(const ColVector & v);
	~ColVector();

	void set(const double v, const int pos);
	double get(const int pos) const;
	int size() const;

	void transpose(RowVector & result) const;

	void neg(ColVector & result) const;
	void neg_assign();

	void mul(ColVector & result, const Matrix & m) const;
	void mul_assign(const Matrix & m);

	ColVector & operator += (const ColVector & v);
	ColVector & operator -= (const ColVector & v);
	ColVector operator + (const ColVector & v) const;
	ColVector operator - (const ColVector & v) const;

	ColVector & operator = (const ColVector & v);

	friend class RowVector;
};

#endif /* MATRIX_H_ */
