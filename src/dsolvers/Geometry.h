/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "flowstar_include.h"
#include "Constraints.h"

class Polyhedron
{
public:
	list<LinearConstraint> constraints;
public:
	Polyhedron();
	Polyhedron(const list<LinearConstraint> & cs);
	Polyhedron(const Polyhedron & P);
	Polyhedron(const Matrix & A, const ColVector & b);
	Polyhedron(const vector<vector<Interval> > & A, const vector<Interval> & B);
	~Polyhedron();

	Interval rho(const vector<Interval> & l) const;
	void tightenConstraints();
	bool empty() const;
	void get(vector<vector<Interval> > & A, vector<Interval> & B) const;

	void dump(FILE *fp, vector<string> const & varNames) const;

	Polyhedron & operator = (const Polyhedron & P);
};

class Parallelotope						// in their constraint-based representations.
{
public:
	Matrix paraTemplate;					// only half of the facet normals are kept
	ColVector b;

	Parallelotope(const Matrix & template_input, const ColVector & b_input);
	Parallelotope(const Parallelotope & P);
	~Parallelotope();

	void center(ColVector & c) const;		// Compute the center point of the parallelotope.
	void dump(FILE *fp) const;

	void toTaylorModel(TaylorModelVec & result) const;		// converse a parallelotope to a Taylor model, the domain is normalized.

	Parallelotope & operator = (const Parallelotope & P);
};

#endif /* GEOMETRY_H_ */
