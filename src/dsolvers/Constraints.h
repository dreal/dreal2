/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#ifndef CONSTRAINTS_H_
#define CONSTRAINTS_H_

#include "Matrix.h"
#include "TaylorModel.h"

// Define the classes of linear and non-linear constraints

class LinearConstraint		// A x <= B not used
{
public:
	vector<Interval> A;
	Interval B;
public:
	LinearConstraint();
	LinearConstraint(const vector<Interval> & A_input, const Interval & B_input);
	LinearConstraint(const LinearConstraint & lc);
	~LinearConstraint();

	void dump(FILE *fp, const vector<string> & stateVarNames) const;

	LinearConstraint & operator = (const LinearConstraint & lc);
};

class PolynomialConstraint	// p(x) <= b
{
public:
	Polynomial p;
	HornerForm hf;		// a HornerForm of p
	Interval B;
public:
	PolynomialConstraint();
	PolynomialConstraint(const Polynomial & p_input, const Interval & B_input);
	PolynomialConstraint(const PolynomialConstraint & pc);
	~PolynomialConstraint();

	void dump(FILE *fp, const vector<string> & stateVarNames) const;

	PolynomialConstraint & operator = (const PolynomialConstraint & pc);
};

#endif /* CONSTRAINTS_H_ */
