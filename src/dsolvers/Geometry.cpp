/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#include "Geometry.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

// class Polyhedron

Polyhedron::Polyhedron()
{
}

Polyhedron::Polyhedron(const list<LinearConstraint> & cs):constraints(cs)
{
}

Polyhedron::Polyhedron(const Polyhedron & P):constraints(P.constraints)
{
}

Polyhedron::Polyhedron(const Matrix & A, const ColVector & b)
{
	int rows = A.rows();
	int cols = A.cols();

	for(int i=0; i<rows; ++i)
	{
		vector<Interval> row;

		for(int j=0; j<cols; ++j)
		{
			Interval intTemp(A.get(i,j));
			row.push_back(intTemp);
		}

		Interval B(b.get(i));
		LinearConstraint lc(row, B);
		constraints.push_back(lc);
	}
}

Polyhedron::Polyhedron(const vector<vector<Interval> > & A, const vector<Interval> & B)
{
	for(int i=0; i<A.size(); ++i)
	{
		LinearConstraint lc(A[i], B[i]);
		constraints.push_back(lc);
	}
}

Polyhedron::~Polyhedron()
{
	constraints.clear();
}

Interval Polyhedron::rho(const vector<Interval> & l) const
{
	int d = l.size();
	int n = constraints.size();
	int size = n*d;

	int *rowInd = new int[ 1 + size ];
	int *colInd = new int[ 1 + size ];
	double *coes = new double [ 1 + size ];

	glp_term_out(GLP_OFF);

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);

	glp_add_rows(lp, n);
	list<LinearConstraint>::const_iterator iter = constraints.begin();
	for(int i=1; i<=n; ++i, ++iter)
		glp_set_row_bnds(lp, i, GLP_UP, 0.0, iter->B.midpoint());

	glp_add_cols(lp, d);
	for(int i=1; i<=d; ++i)
	{
		glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);
		glp_set_obj_coef(lp, i, l[i-1].midpoint());
	}

	iter = constraints.begin();
	for(int i=1; i<=n; ++i, ++iter)
	{
		for(int j=1; j<=d; ++j)
		{
			int pos = j + (i-1)*d;
			rowInd[pos] = i;
			colInd[pos] = j;
			coes[pos] = iter->A[j-1].midpoint();
		}
	}

	glp_load_matrix(lp, size, rowInd, colInd, coes);
	glp_simplex(lp, NULL);
	double result = glp_get_obj_val(lp);
	int status = glp_get_status(lp);

	if(status == GLP_INFEAS || status == GLP_NOFEAS)
		result = INVALID;

	glp_delete_prob(lp);
	delete rowInd;
	delete colInd;
	delete coes;

	Interval intTemp(result);
	return intTemp;
}

void Polyhedron::tightenConstraints()
{
	list<LinearConstraint>::iterator lcIter = constraints.begin();

	for(; lcIter!=constraints.end(); ++lcIter)
	{
		Interval I = rho(lcIter->A);

		if(I < lcIter->B)
		{
			lcIter->B = I;
		}
	}
}

bool Polyhedron::empty() const
{
	if(constraints.size() == 0)
		return false;
	else
	{
		int d = constraints.begin()->A.size();

		vector<Interval> l;
		Interval intZero, intOne(1);
		for(int i=0; i<d; ++i)
		{
			l.push_back(intZero);
		}

		l[0] = intOne;

		Interval result = this->rho(l);

		if(result.inf() <= INVALID + THRESHOLD_HIGH)
			return true;
		else
			return false;
	}
}

void Polyhedron::get(vector<vector<Interval> > & A, vector<Interval> & B) const
{
	A.clear();
	B.clear();

	list<LinearConstraint>::const_iterator iter = constraints.begin();

	for(; iter != constraints.end(); ++iter)
	{
		A.push_back(iter->A);
		B.push_back(iter->B);
	}
}

void Polyhedron::dump(FILE *fp, vector<string> const & varNames) const
{
	list<LinearConstraint>::const_iterator iter = constraints.begin();

	for(; iter != constraints.end(); ++iter)
	{
		iter->dump(fp, varNames);
	}
	fprintf(fp, "\n");
}

Polyhedron & Polyhedron::operator = (const Polyhedron & P)
{
	if(this == &P)
		return *this;

	constraints = P.constraints;
	return *this;
}































// class Parallelotope

Parallelotope::Parallelotope(const Matrix & template_input, const ColVector & b_input):paraTemplate(template_input), b(b_input)
{
}

Parallelotope::Parallelotope(const Parallelotope & P): paraTemplate(P.paraTemplate), b(P.b)
{
}

Parallelotope::~Parallelotope()
{
}

void Parallelotope::center(ColVector & c) const
{
	int d = paraTemplate.cols();

	gsl_vector *r = gsl_vector_alloc(d);
	for(int i=0; i<d; ++i)
		gsl_vector_set( r, i, (b.get(i) - b.get(i+d))/2);

	// We use GSL to solve the linear equations B x = r.

	gsl_matrix *B = gsl_matrix_alloc(d, d);

	for(int i=0; i<d; ++i)
	{
		for(int j=0; j<d; ++j)
		{
			gsl_matrix_set(B, i, j, paraTemplate.get(i,j));
		}
	}

	gsl_vector *x = gsl_vector_alloc(d);

	gsl_linalg_HH_solve(B, r, x);

	for(int i=0; i<d; ++i)
	{
		c.set( gsl_vector_get(x,i), i);
	}

	gsl_vector_free(r);
	gsl_matrix_free(B);
	gsl_vector_free(x);
}

void Parallelotope::dump(FILE *fp) const
{
	int rows = paraTemplate.rows();
	int cols = rows;
	int rangeDim = rows;

	for(int i=0; i<rows; ++i)
	{
		fprintf(fp, "[ ");
		for(int j=0; j<cols-1; ++j)
		{
			fprintf(fp, "%lf,\t", paraTemplate.get(i,j));
		}
		fprintf(fp, "%lf ]\t<=\t%lf\n", paraTemplate.get(i,cols-1), b.get(i));
	}

	for(int i=0; i<rows; ++i)
	{
		fprintf(fp, "[ ");
		for(int j=0; j<cols-1; ++j)
		{
			fprintf(fp, "%lf,\t", -paraTemplate.get(i,j));
		}
		fprintf(fp, "%lf ]\t<=\t%lf\n", -paraTemplate.get(i,cols-1), b.get(i+rangeDim));
	}
}

void Parallelotope::toTaylorModel(TaylorModelVec & result) const
{
	int rangeDim = paraTemplate.rows();
	int domainDim = rangeDim + 1;

	// 1: we converse the center point to a Taylor model

	ColVector colVecCenter(rangeDim);
	center(colVecCenter);

	vector<Interval> coefficients;
	for(int i=0; i<rangeDim; ++i)
	{
		Interval I(colVecCenter.get(i));
		coefficients.push_back(I);
	}

	TaylorModelVec tmvCenter(coefficients, domainDim);

	// 2: we center the parallelotope at 0
	ColVector colVecDiff(rangeDim);
	colVecCenter.mul(colVecDiff, paraTemplate);

	// since a parallelotope is symmetric, we only need to consider half of the intercepts
	ColVector new_b(rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		new_b.set( b.get(i) - colVecDiff.get(i), i);
	}

	// 3: compute the generators.
	Matrix generators(rangeDim, rangeDim);
	vector<int> zeroRows;	// the row indices for zero intercepts

	for(int i=0; i<rangeDim; ++i)
	{
		if(new_b.get(i) <= THRESHOLD_LOW && new_b.get(i) >= -THRESHOLD_LOW)	// zero
		{
			zeroRows.push_back(i);

			for(int j=0; j<rangeDim; ++j)
			{
				generators.set( paraTemplate.get(i,j), i, j);
			}
		}
		else
		{
			for(int j=0; j<rangeDim; ++j)
			{
				generators.set( paraTemplate.get(i,j) / new_b.get(i), i, j);
			}
		}
	}

	generators.inverse_assign();

	Matrix tmv_coefficients(rangeDim, domainDim);

	for(int j=0, k=0; j<rangeDim; ++j)
	{
		if(k < zeroRows.size() && j == zeroRows[k])	// neglect the zero length generators
		{
			++k;
		}
		else
		{
			for(int i=0; i<rangeDim; ++i)
			{
				tmv_coefficients.set( generators.get(i,j), i, j+1);
			}
		}
	}

	TaylorModelVec tmvParallelotope(tmv_coefficients);
	tmvParallelotope.add(result, tmvCenter);
}

Parallelotope & Parallelotope::operator = (const Parallelotope & P)
{
	if(this == &P)
		return *this;

	paraTemplate = P.paraTemplate;
	b = P.b;
	return *this;
}

