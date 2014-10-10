/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#ifndef MONOMIAL_H_
#define MONOMIAL_H_

#include "Interval.h"

extern double cutoff_threshold;

class Monomial
{
private:
	Interval coefficient;		// the coefficient of the monomial
	vector<int> degrees;		// the degrees of the variables, e.g., [2,0,4] is the notation for x1^2 x3^4
	int d;			        	// the degree of the monomial, it is the sum of the values in degrees.

public:
	Monomial();													// empty monomial.
	Monomial(const Interval & I, const vector<int> & degs);
	Monomial(const Monomial & monomial);
	Monomial(const Interval & I, const int numVars);			// a constant
	~Monomial();

	int degree() const;											// degree of the monomial
	int dimension() const;										// dimension of the monomial

	void intEval(Interval & result, const vector<Interval> & domain) const;	// interval evaluation of the monomial

	// interval evaluation of the monomial, we assume that the domain is normalized to [0,s] x [-1,1]^(d-1)
	void intEvalNormal(Interval & result, const vector<Interval> & step_exp_table) const;
	void inv(Monomial & result) const;							// additive inverse

	Monomial & operator = (const Monomial & monomial);
	Monomial & operator += (const Monomial & monomial);			// we assume the two monomials can be added up
	Monomial & operator *= (const Monomial & monomial);
	const Monomial operator + (const Monomial & monomial) const;
	const Monomial operator * (const Monomial & monomial) const;

	bool isLinear(int & index) const;					// Check if the degree of the monomial is 1. If so then return the index of the variable of degree 1.

	void dump_interval(FILE *fp, const vector<string> & varNames) const;	// coefficients are dumped as intervals
	void dump_constant(FILE *fp, const vector<string> & varNames) const;	// coefficients are dumped as constants

	void toString(string & result, const vector<string> & varNames) const;

	bool classInvariantOK() const;

	friend bool operator < (const Monomial & a, const Monomial & b);	// Define a partial order over the monomials
	friend bool operator == (const Monomial & a, const Monomial & b);

	bool cutoff(Monomial & monoRem);
	bool cutoff();

	friend class Polynomial;
	friend class TaylorModel;
	friend class TaylorModelVec;
	friend class Flowpipe;
};


#endif
