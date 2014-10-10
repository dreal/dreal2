/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#include "Polynomial.h"
#include "TaylorModel.h"

vector<Interval> factorial_rec;
vector<Interval> power_4;
vector<Interval> double_factorial;

RangeTree::RangeTree()
{
}

RangeTree::RangeTree(const list<Interval> & ranges_input, const list<RangeTree *> & children_input)
{
	ranges = ranges_input;
	children = children_input;
}

RangeTree::RangeTree(const RangeTree & tree)
{
	ranges = tree.ranges;
	children = tree.children;
}

RangeTree::~RangeTree()
{
	list<RangeTree *>::iterator iter = children.begin();

	for(; iter!=children.end(); ++iter)
	{
		delete *iter;
	}

	ranges.clear();
	children.clear();
}

RangeTree & RangeTree::operator = (const RangeTree & tree)
{
	if(this == &tree)
		return *this;

	ranges = tree.ranges;
	children = tree.children;

	return *this;
}





























// class HornerForm

HornerForm::HornerForm()
{
}

HornerForm::HornerForm(const Interval & I):constant(I)
{
}

HornerForm::HornerForm(const Interval & I, const vector<HornerForm> & hfs):constant(I), hornerForms(hfs)
{
}

HornerForm::HornerForm(const HornerForm & hf):constant(hf.constant), hornerForms(hf.hornerForms)
{
}

HornerForm::~HornerForm()
{
	hornerForms.clear();
}

void HornerForm::clear()
{
	constant.set(0,0);
	hornerForms.clear();
}

void HornerForm::intEval(Interval & result, const vector<Interval> & domain) const
{
	result = constant;

	for(int i=0; i<hornerForms.size(); ++i)
	{
		Interval intHF;
		hornerForms[i].intEval(intHF, domain);
		intHF *= domain[i];
		result += intHF;
	}
}

void HornerForm::insert(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain) const
{
	Interval intZero;
	int numVars = domain.size();

	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert(tmTemp, vars, varsPolyRange, domain);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.remainder *= domain[0];
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert(tmTemp, vars, varsPolyRange, domain);	// recursive call
			tmTemp.mul_insert_assign(vars.tms[i-1], varsPolyRange[i-1], domain);
			result.add_assign(tmTemp);
		}
	}
}

void HornerForm::insert_normal(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars) const
{
	Interval intZero;

	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert_normal(tmTemp, vars, varsPolyRange, step_exp_table, numVars);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.remainder *= step_exp_table[1];
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_normal(tmTemp, vars, varsPolyRange, step_exp_table, numVars);	// recursive call
			tmTemp.mul_insert_normal_assign(vars.tms[i-1], varsPolyRange[i-1], step_exp_table);
			result.add_assign(tmTemp);
		}
	}
}

void HornerForm::insert_ctrunc(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const int order) const
{
	Interval intZero;
	int numVars = domain.size();

	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert_ctrunc(tmTemp, vars, varsPolyRange, domain, order);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.remainder *= domain[0];

		tmTemp.ctrunc(domain, order);
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_ctrunc(tmTemp, vars, varsPolyRange, domain, order);	// recursive call
			tmTemp.mul_insert_ctrunc_assign(vars.tms[i-1], varsPolyRange[i-1], domain, order);
			result.add_assign(tmTemp);
		}
	}
}

void HornerForm::insert_no_remainder(TaylorModel & result, const TaylorModelVec & vars, const int numVars, const int order) const
{
	Interval intZero;
	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert_no_remainder(tmTemp, vars, numVars, order);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.nctrunc(order);
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_no_remainder(tmTemp, vars, numVars, order);	// recursive call
			tmTemp.mul_no_remainder_assign(vars.tms[i-1], order);
			result.add_assign(tmTemp);
		}
	}
}

void HornerForm::insert_no_remainder_no_cutoff(TaylorModel & result, const TaylorModelVec & vars, const int numVars, const int order) const
{
	Interval intZero;
	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert_no_remainder(tmTemp, vars, numVars, order);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.nctrunc(order);
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_no_remainder_no_cutoff(tmTemp, vars, numVars, order);	// recursive call
			tmTemp.mul_no_remainder_no_cutoff_assign(vars.tms[i-1], order);
			result.add_assign(tmTemp);
		}
	}
}

void HornerForm::insert_ctrunc_normal(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const int order) const
{
	Interval intZero;
	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert_ctrunc_normal(tmTemp, vars, varsPolyRange, step_exp_table, numVars, order);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.remainder *= step_exp_table[1];

		tmTemp.ctrunc_normal(step_exp_table, order);
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_ctrunc_normal(tmTemp, vars, varsPolyRange, step_exp_table, numVars, order);	// recursive call

			tmTemp.mul_insert_ctrunc_normal_assign(vars.tms[i-1], varsPolyRange[i-1], step_exp_table, order);

			result.add_assign(tmTemp);
		}
	}
}

void HornerForm::insert_ctrunc_normal(TaylorModel & result, RangeTree * & tree, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const int order) const
{
	Interval intZero;
	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	RangeTree *pnode = new RangeTree;

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		RangeTree *child;

		hornerForms[0].insert_ctrunc_normal(tmTemp, child, vars, varsPolyRange, step_exp_table, numVars, order);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.remainder *= step_exp_table[1];

		Interval intTrunc;
		tmTemp.expansion.ctrunc_normal(intTrunc, step_exp_table, order);
		tmTemp.remainder += intTrunc;

		pnode->ranges.push_back(intTrunc);
		pnode->children.push_back(child);

		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			TaylorModel tmTemp;
			RangeTree *child;

			hornerForms[i].insert_ctrunc_normal(tmTemp, child, vars, varsPolyRange, step_exp_table, numVars, order);	// recursive call

			Interval tm1, intTrunc2;
			tmTemp.mul_insert_ctrunc_normal_assign(tm1, intTrunc2, vars.tms[i-1], varsPolyRange[i-1], step_exp_table, order); 	// here coefficient_range = tm1

			pnode->ranges.push_back(tm1);
			pnode->ranges.push_back(varsPolyRange[i-1]);
			pnode->ranges.push_back(intTrunc2);
			pnode->children.push_back(child);

			result.add_assign(tmTemp);
		}
	}

	tree = pnode;
}

void HornerForm::insert_only_remainder(Interval & result, RangeTree *tree, const TaylorModelVec & vars, const Interval & timeStep) const
{
	Interval intZero;

	result = intZero;
	list<Interval>::const_iterator iter = tree->ranges.begin();
	list<RangeTree *>::const_iterator child = tree->children.begin();

	if(hornerForms.size() > 0)						// the first variable is t
	{
		Interval intTemp;
		hornerForms[0].insert_only_remainder(intTemp, *child, vars, timeStep);
		intTemp *= timeStep;

		intTemp += (*iter);
		result += intTemp;

		++iter;
		++child;

		for(int i=1; i<hornerForms.size(); ++i,++child)
		{
			Interval intTemp2;
			hornerForms[i].insert_only_remainder(intTemp2, *child, vars, timeStep);

			Interval newRemainder = (*iter) * vars.tms[i-1].remainder;
			++iter;
			newRemainder += (*iter) * intTemp2;
			newRemainder += vars.tms[i-1].remainder * intTemp2;
			++iter;
			newRemainder += (*iter);

			result += newRemainder;
			++iter;
		}
	}
}

void HornerForm::dump(FILE *fp, const vector<string> & varNames) const
{
	int numVars = hornerForms.size();

	Interval intZero;
	bool bPlus = false;

	fprintf(fp, " ( ");
	if(!constant.subseteq(intZero))
	{
		bPlus = true;
		constant.dump(fp);
	}

	if(numVars == 0)
	{
		fprintf(fp, " ) ");
		return;
	}

	for(int i=0; i<numVars; ++i)
	{
		if(hornerForms[i].hornerForms.size() != 0 || !hornerForms[i].constant.subseteq(intZero))
		{
			if(bPlus)		// only used to print the "+" symbol
				fprintf(fp, " + ");
			else
				bPlus = true;

			hornerForms[i].dump(fp, varNames);
			fprintf(fp, "* %s", varNames[i].c_str());
		}
	}

	fprintf(fp, " ) ");
}

HornerForm & HornerForm::operator = (const HornerForm & hf)
{
	if(this == &hf)
		return *this;

	constant = hf.constant;
	hornerForms = hf.hornerForms;
	return *this;
}



































// class Polynomial

Polynomial::Polynomial()
{
}

Polynomial::Polynomial(const Interval & constant, const int numVars)
{
	Interval intZero;

	if(!constant.subseteq(intZero))
	{
		Monomial monomial(constant, numVars);
		monomials.push_back(monomial);
	}
}

Polynomial::Polynomial(const RowVector & coefficients)
{
	int numVars = coefficients.size();

	for(int i=0; i<numVars; ++i)
	{
		double dTemp = coefficients.get(i);
		if(dTemp <= THRESHOLD_LOW && dTemp >= -THRESHOLD_LOW)		// dTemp is zero
			continue;

		Interval intTemp(dTemp);
		Monomial monoTemp(intTemp, numVars);
		monoTemp.degrees[i] = 1;
		monoTemp.d = 1;
		monomials.push_back(monoTemp);
	}

	reorder();
}

Polynomial::Polynomial(const vector<Interval> & coefficients)
{
	int numVars = coefficients.size();
	Interval intZero;

	for(int i=0; i<numVars; ++i)
	{
		if(coefficients[i].subseteq(intZero))		// the coefficient is zero
			continue;

		Monomial monoTemp(coefficients[i], numVars);
		monoTemp.degrees[i] = 1;
		monoTemp.d = 1;
		monomials.push_back(monoTemp);
	}

	reorder();
}

Polynomial::Polynomial(const Monomial & monomial)
{
	monomials.push_back(monomial);
}

Polynomial::Polynomial(const list<Monomial> & monos):monomials(monos)
{
	reorder();
}

Polynomial::Polynomial(const Polynomial & polynomial):monomials(polynomial.monomials)
{
}

Polynomial::~Polynomial()
{
	monomials.clear();
}

void Polynomial::reorder()
{
	monomials.sort();
}

void Polynomial::clear()
{
	monomials.clear();
}

void Polynomial::dump_interval(FILE *fp, const vector<string> & varNames) const
{
	if(monomials.size() == 0)
	{
		fprintf(fp, "[0,0]");
		return;
	}

	list<Monomial>::const_iterator iter, iter_last;
	iter_last = monomials.end();
	--iter_last;

	for(iter = monomials.begin(); iter != iter_last; ++iter)
	{
		iter->dump_interval(fp, varNames);
		fprintf(fp, " + ");
	}

	monomials.back().dump_interval(fp, varNames);
}

void Polynomial::dump_constant(FILE *fp, const vector<string> & varNames) const
{
	if(monomials.size() == 0)
	{
		fprintf(fp, "[0,0]");
		return;
	}

	list<Monomial>::const_iterator iter, iter_last;
	iter_last = monomials.end();
	--iter_last;

	for(iter = monomials.begin(); iter != iter_last; ++iter)
	{
		iter->dump_constant(fp, varNames);
		fprintf(fp, " + ");
	}

	monomials.back().dump_constant(fp, varNames);
}

void Polynomial::constant(Interval & result) const
{
	Interval intZero;

	if(monomials.size() > 0 && (monomials.begin())->d == 0)
	{
		result = (monomials.begin())->coefficient;
	}
	else
	{
		result = intZero;
	}
}

void Polynomial::intEval(Interval & result, const vector<Interval> & domain) const
{
	HornerForm hf;
	toHornerForm(hf);
	hf.intEval(result, domain);
}

void Polynomial::intEvalNormal(Interval & result, const vector<Interval> & step_exp_table) const
{
	Interval intZero;
	result = intZero;

	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		Interval intTemp;
		iter->intEvalNormal(intTemp, step_exp_table);

		result += intTemp;
	}
}

void Polynomial::inv(Polynomial & result) const
{
	result = *this;
	result.inv_assign();
}

void Polynomial::inv_assign()
{
	list<Monomial>::iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		iter->coefficient.inv_assign();
	}
}

void Polynomial::add_assign(const Monomial & monomial)
{
	bool bAdded = false;

	list<Monomial>::iterator iter;

	Interval intZero;

	if(monomial.coefficient.subseteq(intZero))
	{
		return;
	}

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		if(monomial < *iter)
		{
			monomials.insert(iter, monomial);
			bAdded = true;
			break;
		}
		else if(monomial == *iter)
		{
			(*iter) += monomial;
/*
			if(iter->coefficient.subseteq(intZero))
			{
				monomials.erase(iter);
			}
*/
			bAdded = true;
			break;
		}
	}

	if(!bAdded)
	{
		monomials.push_back(monomial);
	}
}

void Polynomial::sub_assign(const Monomial & monomial)
{
	Monomial monoTemp;
	monomial.inv(monoTemp);
	add_assign(monoTemp);
}

void Polynomial::mul_assign(const Monomial & monomial)
{
	Interval intZero;

	if(monomial.coefficient.subseteq(intZero))	// the monomial is zero
	{
		clear();
	}
	else
	{
		list<Monomial>::iterator iter;
		for(iter = monomials.begin(); iter != monomials.end(); )
		{
			(*iter) *= monomial;
/*
			if(iter->coefficient.subseteq(intZero))
			{
				iter = monomials.erase(iter);
			}
			else
			{
*/
				++iter;
//			}
		}
	}
}

void Polynomial::mul_assign(const Interval & I)
{
	Interval intZero;

	if(I.subseteq(intZero))	// the interval is zero
	{
		clear();
	}
	else
	{
		list<Monomial>::iterator iter;
		for(iter = monomials.begin(); iter != monomials.end(); )
		{
			iter->coefficient *= I;
			if(iter->coefficient.subseteq(intZero))
			{
				iter = monomials.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
}

void Polynomial::div_assign(const Interval & I)
{
//	Interval intZero;

	list<Monomial>::iterator iter;
	for(iter = monomials.begin(); iter != monomials.end(); )
	{
		iter->coefficient /= I;
/*
		if(iter->coefficient.subseteq(intZero))
		{
			iter = monomials.erase(iter);
		}
		else
		{
*/
			++iter;
//		}
	}
}

void Polynomial::mul(Polynomial & result, const Interval & I) const
{
	result = *this;
	result.mul_assign(I);
}

void Polynomial::div(Polynomial & result, const Interval & I) const
{
	result = *this;
	result.div_assign(I);
}

void Polynomial::mul_assign(const int varIndex, const int degree)
{
	list<Monomial>::iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		iter->degrees[varIndex] += degree;
		iter->d += degree;
	}
}

void Polynomial::mul(Polynomial result, const int varIndex, const int degree) const
{
	result = *this;
	result.mul_assign(varIndex, degree);
}

Polynomial & Polynomial::operator = (const Polynomial & polynomial)
{
	if(this == &polynomial)
		return *this;

	monomials = polynomial.monomials;
	return *this;
}

Polynomial & Polynomial::operator += (const Polynomial & polynomial)
{
	Polynomial result;

//	Interval intZero;

	list<Monomial>::const_iterator iterA;	// polynomial A
	list<Monomial>::const_iterator iterB;	// polynomial B

	for(iterA = monomials.begin(), iterB = polynomial.monomials.begin(); ; )
	{
		if(iterA == monomials.end() || iterB == polynomial.monomials.end())
			break;

		if((*iterA) < (*iterB))
	    {
			result.monomials.push_back(*iterA);
			++iterA;
	    }
		else if((*iterB) < (*iterA))
	    {
			result.monomials.push_back(*iterB);
			++iterB;
	    }
		else
		{
			Interval intTemp;
			intTemp = iterA->coefficient + iterB->coefficient;

//			if(!intTemp.subseteq(intZero))
//			{
				Monomial monoTemp(*iterA);
				monoTemp.coefficient = intTemp;
				result.monomials.push_back(monoTemp);
//			}

			++iterA;
			++iterB;
		}
	}

	if(iterA == monomials.end() && iterB != polynomial.monomials.end())
	{
		for(; iterB != polynomial.monomials.end(); ++iterB)
			result.monomials.push_back(*iterB);
	}
	else if(iterA != monomials.end() && iterB == polynomial.monomials.end())
	{
		for(; iterA != monomials.end(); ++iterA)
			result.monomials.push_back(*iterA);
	}

	*this = result;
	return *this;
}

Polynomial & Polynomial::operator -= (const Polynomial & polynomial)
{
	Polynomial polyTemp = polynomial;
	polyTemp.inv_assign();
	*this += polyTemp;

	return *this;
}

Polynomial & Polynomial::operator *= (const Polynomial & polynomial)
{
	Polynomial result;

	if((monomials.size() == 0) || (polynomial.monomials.size() == 0))
	{
		this->clear();
		return *this;
	}

	list<Monomial>::const_iterator iterB;	// polynomial B

	for(iterB = polynomial.monomials.begin(); iterB != polynomial.monomials.end(); ++iterB)
	{
		Polynomial polyTemp = *this;
		polyTemp.mul_assign(*iterB);
		result += polyTemp;
	}

	*this = result;
	return *this;
}

const Polynomial Polynomial::operator + (const Polynomial & polynomial) const
{
	Polynomial result = *this;
	result += polynomial;
	return result;
}

const Polynomial Polynomial::operator - (const Polynomial & polynomial) const
{
	Polynomial result = *this;
	result -= polynomial;
	return result;
}

const Polynomial Polynomial::operator * (const Polynomial & polynomial) const
{
	Polynomial result = *this;
	result *= polynomial;
	return result;
}

void Polynomial::ctrunc(Interval & remainder, const vector<Interval> & domain, const int order)
{
	Polynomial polyTemp;
	Monomial monoTemp;

	for(; monomials.size() > 0;)
	{
		monoTemp = monomials.back();

		if(monoTemp.d > order)
		{
			polyTemp.monomials.insert(polyTemp.monomials.begin(), monoTemp);
			monomials.pop_back();
		}
		else
		{
			break;
		}
	}

	polyTemp.intEval(remainder, domain);
}

void Polynomial::nctrunc(const int order)
{
	Monomial monoTemp;

	for(; monomials.size() > 0;)
	{
		monoTemp = monomials.back();

		if(monoTemp.d > order)
		{
			monomials.pop_back();
		}
		else
		{
			break;
		}
	}
}

void Polynomial::ctrunc_normal(Interval & remainder, const vector<Interval> & step_exp_table, const int order)
{
	Polynomial polyTemp;
	Monomial monoTemp;

	for(; monomials.size() > 0;)
	{
		monoTemp = monomials.back();

		if(monoTemp.d > order)
		{
			polyTemp.monomials.insert(polyTemp.monomials.begin(), monoTemp);
			monomials.pop_back();
		}
		else
		{
			break;
		}
	}

	polyTemp.intEvalNormal(remainder, step_exp_table);
}

void Polynomial::linearCoefficients(vector<Interval> & result) const
{
	// initially, the result should be filled with 0

	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		int i;

		if(iter->d > 1)
			break;

		if(iter->isLinear(i))
		{
			result[i] = iter->coefficient;
		}
	}
}

void Polynomial::linearCoefficients(RowVector & result) const
{
	// initially, the result should be filled with 0

	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		int i;

		if(iter->d > 1)
			break;

		if(iter->isLinear(i))
		{
			result.set(iter->coefficient.sup(), i);
		}
	}
}

void Polynomial::constraintCoefficients(RowVector & result) const
{
	// initially, the result should be filled with 0

	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		int i;

		if(iter->d > 1)
			break;

		if(iter->isLinear(i))
		{
			if(i > 0)
			{
				result.set(iter->coefficient.sup(), i-1);
			}
		}
	}
}

void Polynomial::constraintCoefficients(vector<Interval> & result) const
{
	// initially, the result should be filled with 0

	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		int i;

		if(iter->d > 1)
			break;

		if(iter->isLinear(i))
		{
			if(i > 0)
			{
				result[i-1] = iter->coefficient;
			}
		}
	}
}

void Polynomial::toHornerForm(HornerForm & result) const
{
	result.clear();

	if(monomials.size() == 0)
		return;

	int numVars = (monomials.begin())->degrees.size();

	list<Monomial> lstMono = monomials;
	list<Monomial>::iterator iter = lstMono.begin();

	if(iter->d == 0)
	{
		result.constant = iter->coefficient;
		iter = lstMono.erase(iter);

		if(lstMono.size() == 0)
			return;
	}

	vector<list<Monomial> > vlMono;

	for(int i=0; i<numVars; ++i)
	{
		list<Monomial> lst_ith;

		for(iter = lstMono.begin(); iter != lstMono.end();)
		{
			if(iter->degrees[i] > 0)
			{
				iter->degrees[i] -= 1;
				iter->d -= 1;
				lst_ith.push_back(*iter);
				iter = lstMono.erase(iter);
			}
			else
			{
				++iter;
			}
		}

		vlMono.push_back(lst_ith);
	}

	for(int i=0; i<numVars; ++i)
	{
		Polynomial polyTemp(vlMono[i]);
		HornerForm hf;
		polyTemp.toHornerForm(hf);
		result.hornerForms.push_back(hf);
	}
}

void Polynomial::rmConstant()
{
	if(monomials.size() > 0 && (monomials.begin())->d == 0)
	{
		monomials.erase( monomials.begin() );
	}
}

int Polynomial::degree() const
{
	if(monomials.size() > 0)
	{
		list<Monomial>::const_iterator iter = monomials.end();
		--iter;
		return iter->d;
	}
	else
	{
		return 0;
	}
}

bool Polynomial::isZero() const
{
	if(monomials.size() == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Polynomial::cutoff_normal(Interval & intRem, const vector<Interval> & step_exp_table)
{
	Polynomial polyTemp;

	list<Monomial>::iterator iter;
	for(iter = monomials.begin(); iter != monomials.end(); )
	{
		Monomial monoTemp;
		bool bvalid = iter->cutoff(monoTemp);

		polyTemp.monomials.push_back(monoTemp);

		if(!bvalid)
		{
			iter = monomials.erase(iter);
		}
		else
		{
			++iter;
		}
	}

	polyTemp.intEvalNormal(intRem, step_exp_table);
}

void Polynomial::cutoff(Interval & intRem, const vector<Interval> & domain)
{
	Polynomial polyTemp;

	list<Monomial>::iterator iter;
	for(iter = monomials.begin(); iter != monomials.end(); )
	{
		Monomial monoTemp;
		bool bvalid = iter->cutoff(monoTemp);

		polyTemp.monomials.push_back(monoTemp);

		if(!bvalid)
		{
			iter = monomials.erase(iter);
		}
		else
		{
			++iter;
		}
	}

	polyTemp.intEval(intRem, domain);
}

void Polynomial::cutoff()
{
	list<Monomial>::iterator iter;
	for(iter = monomials.begin(); iter != monomials.end(); )
	{
		bool bvalid = iter->cutoff();

		if(!bvalid)
		{
			iter = monomials.erase(iter);
		}
		else
		{
			++iter;
		}
	}
}

void Polynomial::derivative(Polynomial & result, const int varIndex) const
{
	result = *this;

	list<Monomial>::iterator iter;

	for(iter = result.monomials.begin(); iter != result.monomials.end(); )
	{
		if(iter->degrees[varIndex] > 0)
		{
			double tmp = iter->degrees[varIndex];
			iter->degrees[varIndex] -= 1;
			iter->d -= 1;
			iter->coefficient.mul_assign(tmp);
			++iter;
		}
		else
		{
			iter = result.monomials.erase(iter);
		}
	}
}

void Polynomial::LieDerivative(Polynomial & result, const vector<Polynomial> & f) const
{
	derivative(result, 0);

	int rangeDim = f.size();

	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial P;
		derivative(P, i+1);
		P *= f[i];
		result += P;
	}
}

void Polynomial::sub(Polynomial & result, const Polynomial & P, const int order) const
{
	list<Monomial> monomials1, monomials2;
	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		if(iter->d == order)
		{
			monomials1.push_back(*iter);
		}
	}

	for(iter = P.monomials.begin(); iter != P.monomials.end(); ++iter)
	{
		if(iter->d == order)
		{
			monomials2.push_back(*iter);
		}
	}

	Polynomial P1(monomials1), P2(monomials2);
	result = P1 - P2;
}

void Polynomial::exp_taylor(Polynomial & result, const int numVars, const int order) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();				// F = tm - c

	const_part.exp_assign();	// exp(c)

	if(F.isZero())				// tm = c
	{
		Polynomial polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	Interval I(1);
	Polynomial polyOne(I, numVars);

	// to compute the expression 1 + F + (1/2!)F^2 + ... + (1/k!)F^k,
	// we evaluate its Horner form (...((1/(k-1))((1/k)*F+1)*F + 1) ... + 1)

	result = polyOne;

	for(int i=order; i>0; --i)
	{
		Interval intFactor(1);
		intFactor.div_assign((double)i);

		result.mul_assign(intFactor);

		result *= F;
		result.nctrunc(order);
		result.cutoff();

		result += polyOne;
	}

	result.mul_assign(const_part);
}

void Polynomial::rec_taylor(Polynomial & result, const int numVars, const int order) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();				// F = tm - c

	const_part.rec_assign();	// 1/c

	if(F.isZero())				// tm = c
	{
		Polynomial polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	Interval I(1);
	Polynomial polyOne(I, numVars);
	Polynomial F_c;
	F.mul(F_c, const_part);

	// to compute the expression 1 - F/c + (F/c)^2 - ... + (-1)^k (F/c)^k,
	// we evaluate its Horner form (-1)*(...((-1)*(-F/c + 1)*F/c + 1)...) + 1

	result = polyOne;

	for(int i=order; i>0; --i)
	{
		result.inv_assign();

		result *= F_c;
		result.nctrunc(order);
		result.cutoff();

		result += polyOne;
	}

	result.mul_assign(const_part);
}

void Polynomial::sin_taylor(Polynomial & result, const int numVars, const int order) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();			// F = tm - c

	if(F.isZero())			// tm = c
	{
		const_part.sin_assign();
		Polynomial polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	Interval sinc, cosc, msinc, mcosc;
	sinc = const_part.sin();
	cosc = const_part.cos();
	sinc.inv(msinc);
	cosc.inv(mcosc);

	Polynomial polyTemp(sinc, numVars);
	result = polyTemp;

	int k=1;
	Interval I(1);

	Polynomial polyPowerF(I, numVars);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		I.div_assign((double)i);

		switch(k)
		{
		case 0:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff();

			polyTemp = polyPowerF;
			Interval intTemp = I * sinc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 1:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff();

			polyTemp = polyPowerF;
			Interval intTemp = I * cosc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 2:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff();

			polyTemp = polyPowerF;
			Interval intTemp = I * msinc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 3:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff();

			polyTemp = polyPowerF;
			Interval intTemp = I * mcosc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		}
	}

	result.cutoff();
}

void Polynomial::cos_taylor(Polynomial & result, const int numVars, const int order) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();			// F = tm - c

	if(F.isZero())			// tm = c
	{
		const_part.cos_assign();
		Polynomial polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	Interval sinc, cosc, msinc, mcosc;
	sinc = const_part.sin();
	cosc = const_part.cos();
	sinc.inv(msinc);
	cosc.inv(mcosc);

	Polynomial polyTemp(cosc, numVars);
	result = polyTemp;

	int k=1;
	Interval I(1);

	Polynomial polyPowerF(I, numVars);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		I.div_assign((double)i);

		switch(k)
		{
		case 0:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff();

			polyTemp = polyPowerF;
			Interval intTemp = I * cosc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 1:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff();

			polyTemp = polyPowerF;
			Interval intTemp = I * msinc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 2:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff();

			polyTemp = polyPowerF;
			Interval intTemp = I * mcosc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 3:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff();

			polyTemp = polyPowerF;
			Interval intTemp = I * sinc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		}
	}

	result.cutoff();
}

void Polynomial::log_taylor(Polynomial & result, const int numVars, const int order) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();			// F = tm - c

	Interval C = const_part;

	const_part.log_assign();	// log(c)

	if(F.isZero())			// tm = c
	{
		Polynomial polyLog(const_part, numVars);
		result = polyLog;

		return;
	}

	Polynomial F_c;
	F.div(F_c, C);

	result = F_c;

	Interval I((double)order);
	result.div_assign(I);			// F/c * (1/order)

	for(int i=order; i>=2; --i)
	{
		Interval J(1);
		J.div_assign((double)(i-1));
		Polynomial polyJ(J, numVars);

		result -= polyJ;
		result.inv_assign();

		result *= F_c;
		result.nctrunc(order);
		result.cutoff();
	}

	Polynomial const_part_poly(const_part, numVars);
	result += const_part_poly;
}

void Polynomial::sqrt_taylor(Polynomial & result, const int numVars, const int order) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();			// F = tm - c

	Interval C = const_part;
	const_part.sqrt_assign();	// log(c)

	if(F.isZero())			// tm = c
	{
		Polynomial polySqrt(const_part, numVars);
		result = polySqrt;

		return;
	}

	Polynomial F_2c;
	F.div(F_2c, C);

	Interval intTwo(2);
	F_2c.div_assign(intTwo);	// F/2c

	Interval intOne(1);
	Polynomial polyOne(intOne, numVars);

	result = F_2c;

	Interval K(1), J(1);

	for(int i=order, j=2*order-3; i>=2; --i, j-=2)
	{
		// i
		Interval K((double)i);

		// j = 2*i-3
		Interval J((double)j);

		result.inv_assign();
		result.mul_assign( J / K );

		result += polyOne;
		result *= F_2c;
		result.nctrunc(order);
		result.cutoff();
	}

	result += polyOne;

	result.mul_assign(const_part);
}

void Polynomial::toString(string & result, const vector<string> & varNames) const
{
	string strPoly;

	if(monomials.size() == 0)
	{
		strPoly = "(0)";
		return;
	}

	list<Monomial>::const_iterator iter, iter_last;
	iter_last = monomials.end();
	--iter_last;

	strPoly += '(';

	for(iter = monomials.begin(); iter != iter_last; ++iter)
	{
		string strTemp;
		iter->toString(strTemp, varNames);

		strPoly += strTemp;
		strPoly += ' ';
		strPoly += '+';
		strPoly += ' ';
	}

	string strTemp2;
	monomials.back().toString(strTemp2, varNames);
	strPoly += strTemp2;
	strPoly += ')';

	result = strPoly;
}









































void compute_factorial_rec(const int order)
{
	Interval I(1);

	factorial_rec.push_back(I);

	for(int i=1; i<=order; ++i)
	{
		I.div_assign((double)i);
		factorial_rec.push_back(I);
	}
}

void compute_power_4(const int order)
{
	Interval I(1);

	power_4.push_back(I);

	for(int i=1; i<=order; ++i)
	{
		I.mul_assign(4.0);
		power_4.push_back(I);
	}
}

void compute_double_factorial(const int order)
{
	Interval odd(1), even(1);

	double_factorial.push_back(even);
	double_factorial.push_back(odd);

	for(int i=2; i<=order; ++i)
	{
		if(i%2 == 0)
		{
			even.mul_assign((double)i);
			double_factorial.push_back(even);
		}
		else
		{
			odd.mul_assign((double)i);
			double_factorial.push_back(odd);
		}
	}
}

void computeTaylorExpansion(vector<HornerForm> & result, const vector<Polynomial> & ode, const int order)
{
	int rangeDim = ode.size();

	vector<Polynomial> taylorExpansion;
	vector<Polynomial> LieDeriv_n;

	for(int i=0; i<rangeDim; ++i)
	{
		RowVector row(rangeDim+1);
		row.set(1,i+1);
		Polynomial P(row);
		taylorExpansion.push_back(P);
		LieDeriv_n.push_back(P);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=1; j<=order; ++j)
		{
			Polynomial P1;
			LieDeriv_n[i].LieDerivative(P1, ode);
			LieDeriv_n[i] = P1;

			P1.mul_assign(factorial_rec[j]);
			P1.mul_assign(0,j);

			taylorExpansion[i] += P1;
		}

		taylorExpansion[i].cutoff();
	}

	result.clear();

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		HornerForm hf;
		taylorExpansion[i].toHornerForm(hf);
		result.push_back(hf);
	}
}

void computeTaylorExpansion(vector<HornerForm> & result, const vector<Polynomial> & ode, const vector<int> & orders)
{
	int rangeDim = ode.size();

	vector<Polynomial> taylorExpansion;
	vector<Polynomial> LieDeriv_n;

	for(int i=0; i<rangeDim; ++i)
	{
		RowVector row(rangeDim+1);
		row.set(1,i+1);
		Polynomial P(row);
		taylorExpansion.push_back(P);
		LieDeriv_n.push_back(P);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=1; j<=orders[i]; ++j)
		{
			Polynomial P1;
			LieDeriv_n[i].LieDerivative(P1, ode);
			LieDeriv_n[i] = P1;

			P1.mul_assign(factorial_rec[j]);
			P1.mul_assign(0,j);

			taylorExpansion[i] += P1;
		}

		taylorExpansion[i].cutoff();
	}

	result.clear();

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		HornerForm hf;
		taylorExpansion[i].toHornerForm(hf);
		result.push_back(hf);
	}
}

void computeTaylorExpansion(vector<HornerForm> & resultHF, vector<Polynomial> & resultMF, vector<Polynomial> & highest, const vector<Polynomial> & ode, const int order)
{
	int rangeDim = ode.size();

	vector<Polynomial> taylorExpansion;
	vector<Polynomial> LieDeriv_n;

	for(int i=0; i<rangeDim; ++i)
	{
		RowVector row(rangeDim+1);
		row.set(1,i+1);
		Polynomial P(row);
		taylorExpansion.push_back(P);
		LieDeriv_n.push_back(P);
	}

	highest.clear();

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=1; j<=order; ++j)
		{
			Polynomial P;
			LieDeriv_n[i].LieDerivative(P, ode);
			LieDeriv_n[i] = P;

			if(j == order)
			{
				highest.push_back(P);
			}

			P.mul_assign(factorial_rec[j]);
			P.mul_assign(0,j);

			taylorExpansion[i] += P;
		}

		taylorExpansion[i].cutoff();
	}

	resultMF = taylorExpansion;

	resultHF.clear();
	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		HornerForm hf;
		taylorExpansion[i].toHornerForm(hf);
		resultHF.push_back(hf);
	}
}

void computeTaylorExpansion(vector<HornerForm> & resultHF, vector<Polynomial> & resultMF, vector<Polynomial> & highest, const vector<Polynomial> & ode, const vector<int> & orders)
{
	int rangeDim = ode.size();

	vector<Polynomial> taylorExpansion;
	vector<Polynomial> LieDeriv_n;

	for(int i=0; i<rangeDim; ++i)
	{
		RowVector row(rangeDim+1);
		row.set(1,i+1);
		Polynomial P(row);
		taylorExpansion.push_back(P);
		LieDeriv_n.push_back(P);
	}

	highest.clear();

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=1; j<=orders[i]; ++j)
		{
			Polynomial P;
			LieDeriv_n[i].LieDerivative(P, ode);
			LieDeriv_n[i] = P;

			if(j == orders[i])
			{
				highest.push_back(P);
			}

			P.mul_assign(factorial_rec[j]);
			P.mul_assign(0,j);

			taylorExpansion[i] += P;
		}

		taylorExpansion[i].cutoff();
	}

	resultMF = taylorExpansion;

	resultHF.clear();
	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		HornerForm hf;
		taylorExpansion[i].toHornerForm(hf);
		resultHF.push_back(hf);
	}
}

void increaseExpansionOrder(vector<HornerForm> & resultHF, vector<Polynomial> & resultMF, vector<Polynomial> & highest, const vector<Polynomial> & taylorExpansion, const vector<Polynomial> & ode, const int order)
{
	int rangeDim = ode.size();

	vector<Polynomial> expansion = taylorExpansion;
	vector<Polynomial> LieDeriv_n = highest;

	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial P1;
		LieDeriv_n[i].LieDerivative(P1, ode);

		highest[i] = P1;

		P1.mul_assign(factorial_rec[order+1]);
		P1.mul_assign(0, order+1);

		expansion[i] += P1;

		expansion[i].cutoff();
	}

	resultMF = expansion;

	resultHF.clear();
	for(int i=0; i<expansion.size(); ++i)
	{
		HornerForm hf;
		expansion[i].toHornerForm(hf);
		resultHF.push_back(hf);
	}
}

void increaseExpansionOrder(HornerForm & resultHF, Polynomial & resultMF, Polynomial & highest, const Polynomial & taylorExpansion, const vector<Polynomial> & ode, const int order)
{
	int rangeDim = ode.size();

	Polynomial expansion = taylorExpansion;
	Polynomial LieDeriv_n = highest;

	Polynomial P1;
	LieDeriv_n.LieDerivative(P1, ode);

	highest = P1;

	P1.mul_assign(factorial_rec[order+1]);
	P1.mul_assign(0, order+1);

	expansion += P1;
	expansion.cutoff();

	resultMF = expansion;
	expansion.toHornerForm(resultHF);
}




