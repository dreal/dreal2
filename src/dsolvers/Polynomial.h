/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include "Monomial.h"
#include "Matrix.h"

class Polynomial;
class TaylorModel;
class TaylorModelVec;
class Flowpipe;

extern vector<Interval> factorial_rec;
extern vector<Interval> power_4;
extern vector<Interval> double_factorial;

class RangeTree
{
public:
	list<Interval> ranges;
	list<RangeTree *> children;

	RangeTree();
	RangeTree(const list<Interval> & ranges_input, const list<RangeTree *> & children_input);
	RangeTree(const RangeTree & tree);
	~RangeTree();

	RangeTree & operator = (const RangeTree & tree);
};

class HornerForm							// c + (...)*x1 + (...)*x2 + ... + (...)*xn
{
private:
	Interval constant;						// constant part
	vector<HornerForm> hornerForms;			// other parts
public:
	HornerForm();
	HornerForm(const Interval & I);
	HornerForm(const Interval & I, const vector<HornerForm> & hfs);
	HornerForm(const HornerForm & hf);
	~HornerForm();

	void clear();
	void intEval(Interval & result, const vector<Interval> & domain) const;		// interval evaluation of the Horner form

	// substitute the variables by the given Taylor models
	void insert(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain) const;
	void insert_normal(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars) const;

	// with conservative truncation
	void insert_ctrunc(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const int order) const;
	// with non-conservative truncation
	void insert_no_remainder(TaylorModel & result, const TaylorModelVec & vars, const int numVars, const int order) const;
	void insert_no_remainder_no_cutoff(TaylorModel & result, const TaylorModelVec & vars, const int numVars, const int order) const;

	void insert_ctrunc_normal(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const int order) const;

	void insert_ctrunc_normal(TaylorModel & result, RangeTree * & tree, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const int order) const;	// the first time
	void insert_only_remainder(Interval & result, RangeTree *tree, const TaylorModelVec & vars, const Interval & timeStep) const;	// after the first time

	void dump(FILE *fp, const vector<string> & varNames) const;	// only for tests

	HornerForm & operator = (const HornerForm & hf);

	friend class Polynomial;
};

class Polynomial				// polynomials in monomial form
{
private:
	list<Monomial> monomials;
public:
	Polynomial();														// empty polynomial
	Polynomial(const Interval & constant, const int numVars);			// constant polynomial where dim is the number of the variables
	Polynomial(const RowVector & coefficients);
	Polynomial(const vector<Interval> & coefficients);					// linear polynomial with the given coefficients, the input matrix is a row vector
	Polynomial(const Monomial & monomial);								// polynomial with one monomial
	Polynomial(const list<Monomial> & monos);
	Polynomial(const Polynomial & polynomial);
	~Polynomial();

	void reorder();														// sort the monomials.
	void clear();

	void dump_interval(FILE *fp, const vector<string> & varNames) const;
	void dump_constant(FILE *fp, const vector<string> & varNames) const;

	void constant(Interval & result) const;											// constant part of the polynomial
	void intEval(Interval & result, const vector<Interval> & domain) const;			// interval evaluation of the polynomial
	void intEvalNormal(Interval & result, const vector<Interval> & step_exp_table) const;	// fast evaluation over normalized domain
	void inv(Polynomial & result) const;											// additive inverse
	void inv_assign();
	void add_assign(const Monomial & monomial);										// add a monomial
	void sub_assign(const Monomial & monomial);										// subtract a monomial
	void mul_assign(const Monomial & monomial);										// multiplied by a monomial

	void mul_assign(const Interval & I);											// multiplied by an interval
	void div_assign(const Interval & I);											// divided by an interval
	void mul(Polynomial & result, const Interval & I) const;
	void div(Polynomial & result, const Interval & I) const;

	void mul_assign(const int varIndex, const int degree);							// multiplied by a term x^d
	void mul(Polynomial result, const int varIndex, const int degree) const;

	Polynomial & operator = (const Polynomial & P);
	Polynomial & operator += (const Polynomial & polynomial);
	Polynomial & operator -= (const Polynomial & polynomial);
	Polynomial & operator *= (const Polynomial & polynomial);
	const Polynomial operator + (const Polynomial & polynomial) const;
	const Polynomial operator - (const Polynomial & polynomial) const;
	const Polynomial operator * (const Polynomial & polynomial) const;

	void ctrunc(Interval & remainder, const vector<Interval> & domain, const int order);	// conservative truncation
	void nctrunc(const int order);															// non-conservative truncation
	void ctrunc_normal(Interval & remainder, const vector<Interval> & step_exp_table, const int order);

	void linearCoefficients(vector<Interval> & result) const;								// the coefficients of the linear part, e.g. x+2y+z^3 returns [1,2,0]
	void linearCoefficients(RowVector & result) const;
	void constraintCoefficients(RowVector & result) const;
	void constraintCoefficients(vector<Interval> & result) const;
	void toHornerForm(HornerForm & result) const;											// transform the polynomial into a Horner form

	void rmConstant();				// remove the constant part
	int degree() const;				// degree of the polynomial
	bool isZero() const;

	void cutoff_normal(Interval & intRem, const vector<Interval> & step_exp_table);
	void cutoff(Interval & intRem, const vector<Interval> & domain);
	void cutoff();

	void derivative(Polynomial & result, const int varIndex) const;					// derivative with respect to a variable
	void LieDerivative(Polynomial & result, const vector<Polynomial> & f) const;	// Lie derivative without truncation

	void sub(Polynomial & result, const Polynomial & P, const int order) const;		// compute the subtraction of the monomials with some order

	void exp_taylor(Polynomial & result, const int numVars, const int order) const;
	void rec_taylor(Polynomial & result, const int numVars, const int order) const;
	void sin_taylor(Polynomial & result, const int numVars, const int order) const;
	void cos_taylor(Polynomial & result, const int numVars, const int order) const;
	void log_taylor(Polynomial & result, const int numVars, const int order) const;
	void sqrt_taylor(Polynomial & result, const int numVars, const int order) const;

	void toString(string & result, const vector<string> & varNames) const;	// transform a polynomial to a string

	friend class TaylorModel;
	friend class TaylorModelVec;
	friend class Flowpipe;
	friend class ContinuousSystem;
};

void compute_factorial_rec(const int order);
void compute_power_4(const int order);
void compute_double_factorial(const int order);

void computeTaylorExpansion(vector<HornerForm> & result, const vector<Polynomial> & ode, const int order);
void computeTaylorExpansion(vector<HornerForm> & result, const vector<Polynomial> & ode, const vector<int> & orders);

void computeTaylorExpansion(vector<HornerForm> & resultHF, vector<Polynomial> & resultMF, vector<Polynomial> & highest, const vector<Polynomial> & ode, const int order);
void computeTaylorExpansion(vector<HornerForm> & resultHF, vector<Polynomial> & resultMF, vector<Polynomial> & highest, const vector<Polynomial> & ode, const vector<int> & orders);

void increaseExpansionOrder(vector<HornerForm> & resultHF, vector<Polynomial> & resultMF, vector<Polynomial> & highest, const vector<Polynomial> & taylorExpansion, const vector<Polynomial> & ode, const int order);
void increaseExpansionOrder(HornerForm & resultHF, Polynomial & resultMF, Polynomial & highest, const Polynomial & taylorExpansion, const vector<Polynomial> & ode, const int order);

#endif /* POLYNOMIAL_H_ */
