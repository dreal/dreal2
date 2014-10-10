/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#include "Constraints.h"

// class linearConstraint

LinearConstraint::LinearConstraint()
{
}

LinearConstraint::LinearConstraint(const vector<Interval> & A_input, const Interval & B_input)
{
	A = A_input;
	B = B_input;
}

LinearConstraint::LinearConstraint(const LinearConstraint & lc)
{
	A = lc.A;
	B = lc.B;
}

LinearConstraint::~LinearConstraint()
{
	A.clear();
}

void LinearConstraint::dump(FILE *fp, const vector<string> & stateVarNames) const
{
	int d = A.size();
	for(int i=0; i<d-1; ++i)
	{
		fprintf(fp, "(%lf*%s) + ", A[i].midpoint(), stateVarNames[i].c_str());
	}

	fprintf(fp, "(%lf*%s)", A[d-1].midpoint(), stateVarNames[d-1].c_str());

	fprintf(fp, " <= %lf\n", B.midpoint());
}

LinearConstraint & LinearConstraint::operator = (const LinearConstraint & lc)
{
	if(this == &lc)
		return *this;

	A = lc.A;
	B = lc.B;
	return *this;
}


































// class PolynomialConstraint

PolynomialConstraint::PolynomialConstraint()
{
}

PolynomialConstraint::PolynomialConstraint(const Polynomial & p_input, const Interval & B_input)
{
	p = p_input;
	p.toHornerForm(hf);
	B = B_input;
}

PolynomialConstraint::PolynomialConstraint(const PolynomialConstraint & pc)
{
	p = pc.p;
	hf = pc.hf;
	B = pc.B;
}

PolynomialConstraint::~PolynomialConstraint()
{
}

void PolynomialConstraint::dump(FILE *fp, const vector<string> & stateVarNames) const
{
	string tVar("local_t");
	vector<string> names = stateVarNames;
	names.insert(names.begin(), tVar);

	p.dump_constant(fp, names);
	fprintf(fp, " <= %lf\n", B.midpoint());
}

PolynomialConstraint & PolynomialConstraint::operator = (const PolynomialConstraint & pc)
{
	if(this == &pc)
		return *this;

	p = pc.p;
	hf = pc.hf;
	B = pc.B;

	return *this;
}
