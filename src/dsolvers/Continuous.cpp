/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#include "Continuous.h"

Flowpipe::Flowpipe()
{
}

Flowpipe::Flowpipe(const TaylorModelVec & tmvPre_input, const TaylorModelVec & tmv_input, const vector<Interval> & domain_input):
		tmvPre(tmvPre_input), tmv(tmv_input), domain(domain_input)
{
}

Flowpipe::Flowpipe(const vector<Interval> & box, const Interval & I)
{
	int rangeDim = box.size();
	int domainDim = rangeDim + 1;
	Interval intUnit(-1,1), intZero;

	TaylorModelVec tmvCenter;
	vector<double> scalars;

	domain.push_back(I);		// time interval

	// normalize the domain to [-1,1]^n
	for(int i=0; i<rangeDim; ++i)
	{
		double midpoint = box[i].midpoint();
		Interval intMid(midpoint);
		TaylorModel tmTemp(intMid, domainDim);
		tmvCenter.tms.push_back(tmTemp);

		Interval intTemp = box[i];
		intTemp.sub_assign(midpoint);
		scalars.push_back( intTemp.sup() );
		domain.push_back(intUnit);
	}

	Matrix coefficients_of_tmvPre(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		coefficients_of_tmvPre.set(scalars[i], i, i+1);
	}

	TaylorModelVec tmvTemp(coefficients_of_tmvPre);
	tmvTemp.add(tmvPre, tmvCenter);

	Matrix coefficients_of_tmv(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		coefficients_of_tmv.set(1, i, i+1);
	}

	TaylorModelVec tmvTemp2(coefficients_of_tmv);
	tmv = tmvTemp2;
}

Flowpipe::Flowpipe(const TaylorModelVec & tmv_input, const vector<Interval> & domain_input)
{
	int rangeDim = tmv_input.tms.size();

	tmv = tmv_input;
	domain = domain_input;

	Matrix coefficients_of_tmvPre(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		coefficients_of_tmvPre.set(1, i, i+1);
	}

	TaylorModelVec tmvTemp(coefficients_of_tmvPre);
	tmvPre = tmvTemp;

	normalize();
}

Flowpipe::Flowpipe(const Flowpipe & flowpipe):tmvPre(flowpipe.tmvPre), tmv(flowpipe.tmv), domain(flowpipe.domain)
{
}

Flowpipe::~Flowpipe()
{
	clear();
}

void Flowpipe::clear()
{
	tmvPre.clear();
	tmv.clear();
	domain.clear();
}

void Flowpipe::dump(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames) const
{
	TaylorModelVec tmvTemp;

	composition(tmvTemp);

	// dump the Taylor model
	tmvTemp.dump_interval(fp, stateVarNames, tmVarNames);

	//dump the domain
	for(int i=0; i<domain.size(); ++i)
	{
		fprintf(fp, "%s in ", tmVarNames[i].c_str());
		domain[i].dump(fp);
		fprintf(fp, "\n");
	}
}

void Flowpipe::dump_normal(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames, vector<Interval> & step_exp_table) const
{
	TaylorModelVec tmvTemp;

	composition_normal(tmvTemp, step_exp_table);

	// dump the Taylor model
	tmvTemp.dump_interval(fp, stateVarNames, tmVarNames);

	//dump the domain
	for(int i=0; i<domain.size(); ++i)
	{
		fprintf(fp, "%s in ", tmVarNames[i].c_str());
		domain[i].dump(fp);
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n");
}

void Flowpipe::composition(TaylorModelVec & result) const
{
	vector<int> orders;

	for(int i=0; i<tmv.tms.size(); ++i)
	{
		int d1 = tmv.tms[i].degree();
		int d2 = tmvPre.tms[i].degree();

		if(d1 > d2)
		{
			orders.push_back(d1);
		}
		else
		{
			orders.push_back(d2);
		}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRange(tmvPolyRange, domain);
	tmvPre.insert_ctrunc(result, tmv, tmvPolyRange, domain, orders);
}

void Flowpipe::composition_normal(TaylorModelVec & result, const vector<Interval> & step_exp_table) const
{
	int domainDim = domain.size();

	vector<int> orders;

	for(int i=0; i<tmv.tms.size(); ++i)
	{
		int d1 = tmv.tms[i].degree();
		int d2 = tmvPre.tms[i].degree();

		if(d1 > d2)
		{
			orders.push_back(d1);
		}
		else
		{
			orders.push_back(d2);
		}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_exp_table);
	tmvPre.insert_ctrunc_normal(result, tmv, tmvPolyRange, step_exp_table, domainDim, orders);
}

void Flowpipe::intEval(vector<Interval> & result) const
{
	TaylorModelVec tmvTemp;
	composition(tmvTemp);

	tmvTemp.intEval(result, domain);
}

void Flowpipe::intEvalNormal(vector<Interval> & result, const vector<Interval> & step_exp_table) const
{
	TaylorModelVec tmvTemp;
	composition_normal(tmvTemp, step_exp_table);
	tmvTemp.intEvalNormal(result, step_exp_table);
}

void Flowpipe::normalize()
{
	Interval intZero;

	// we first normalize the Taylor model tmv
	tmv.normalize(domain);

	int rangeDim = tmv.tms.size();

	// compute the center point of tmv
	vector<Interval> intVecCenter;
	tmv.constant(intVecCenter);
	tmv.rmConstant();

	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		tmv.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	vector<Interval> tmvRange;
	tmv.intEval(tmvRange, domain);

	vector<vector<Interval> > coefficients;
	vector<Interval> row;

	for(int i=0; i<rangeDim+1; ++i)
	{
		row.push_back(intZero);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		coefficients.push_back(row);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		Interval intScalor;
		tmvRange[i].mag(intScalor);

		if(intScalor.subseteq(intZero))
		{
			coefficients[i][i+1] = intZero;
		}
		else
		{
			coefficients[i][i+1] = intScalor;
			tmv.tms[i].div_assign(intScalor);
		}
	}

	TaylorModelVec newVars(coefficients);
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp(intVecCenter[i], rangeDim+1);
		newVars.tms[i].add_assign(tmTemp);
	}

	for(int i=0; i<tmvPre.tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tmvPre.tms[i].insert_no_remainder_no_cutoff(tmTemp, newVars, rangeDim+1, tmvPre.tms[i].degree());
		tmvPre.tms[i].expansion = tmTemp.expansion;
	}
}

// for low-degree ODEs
// fixed step sizes and orders

bool Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);

	// the center point of the remainder
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), order);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, order);
		x.tms.push_back(tmTemp);
	}

	x.cutoff();

	bool bfound = true;

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];	// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);

		intDifferences.push_back(intTemp);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return false;
	}
	else
	{
		for(int i=0; i<rangeDim; ++i)
			x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders.push_back(tmvTemp.tms[i].remainder);
		}

		// add the uncertainties and the cutoff intervals onto the result
		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}

bool Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, orders[i]);
		x.tms.push_back(tmTemp);
	}

	x.cutoff();

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return false;
	}
	else
	{
		for(int i=0; i<rangeDim; ++i)
			x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}

// adaptive step sizes and fixed orders

bool Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), order);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, order);
		x.tms.push_back(tmTemp);
	}

	x.cutoff();

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
	}

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return false;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*order);

		for(int i=0; i<rangeDim; ++i)	// update the uncertainties in the step
		{
			step_uncertainties[i] = step_exp_table[1] * uncertainties[i];
		}

		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}

bool Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, orders[i]);
		x.tms.push_back(tmTemp);
	}

	x.cutoff();

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);
	}

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return false;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*globalMaxOrder);

		for(int i=0; i<rangeDim; ++i)	// update the uncertainties in the step
		{
			step_uncertainties[i] = step_exp_table[1] * uncertainties[i];
		}

		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders);

		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}

// adaptive orders and fixed step sizes

bool Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), order);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, order);
		x.tms.push_back(tmTemp);
	}

	x.cutoff();

	bool bfound = true;

	vector<Interval> step_uncertainties;

	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	// add the uncertainties onto the reault
	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	int newOrder = order;

	for(; !bfound;)
	{
		++newOrder;

		if(newOrder > maxOrder)
		{
			return false;
		}

		// increase the approximation orders by 1
		x.Picard_no_remainder_assign(x0, ode, rangeDim+1, newOrder);

		for(int i=0; i<rangeDim; ++i)	// apply the estimation again
		{
			x.tms[i].remainder = estimation[i] + step_uncertainties[i];
		}

		// compute the Picard operation again
		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, newOrder);

		// Update the irreducible part
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial polyTemp;
			polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

			Interval intTemp;
			polyTemp.intEvalNormal(intTemp, step_exp_table);
			intDifferences[i] = intTemp;
		}

		// add the uncertainties onto the result
		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;
		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		// add the uncertainties
		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	order = newOrder;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}

bool Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, orders[i]);
		x.tms.push_back(tmTemp);
	}

	x.cutoff();

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	vector<bool> bIncrease;
	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bIncrease.push_back(true);
			if(bfound)
				bfound = false;
		}
		else
		{
			bIncrease.push_back(false);
		}
	}

	vector<int> newOrders = orders;
	bool bIncreaseOthers = false;
	int numIncrease = 0;

	vector<bool> bIncreased;
	for(int i=0; i<rangeDim; ++i)
	{
		bIncreased.push_back(false);
	}

	for(; !bfound;)
	{
		bool bChanged = false;

		if(bIncreaseOthers)
		{
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(!bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			bIncreaseOthers = false;
		}
		else
		{
			numIncrease = 0;
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					++numIncrease;

					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			if(numIncrease < newOrders.size())
				bIncreaseOthers = true;
		}

		if(!bChanged)
		{
			return false;
		}

		// increase the approximation orders
		x.Picard_no_remainder_assign(x0, ode, rangeDim+1, newOrders, bIncreased);

		for(int i=0; i<rangeDim; ++i)	// apply the estimation again
		{
			x.tms[i].remainder = estimation[i] + step_uncertainties[i];
		}

		// compute the Picard operation again
		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, newOrders);

		for(int i=0; i<rangeDim; ++i)
		{
			// update the irreducible part if necessary
			if(bIncreased[i])
			{
				Polynomial polyTemp;
				polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

				Interval intTemp;
				polyTemp.intEvalNormal(intTemp, step_exp_table);
				intDifferences[i] = intTemp;
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = true;
					if(bfound)
						bfound = false;
				}
				else
				{
					bfound = false;
					break;
				}
			}
			else
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = false;
				}
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			bIncreased[i] = false;
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	orders = newOrders;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}























// for high-degree ODEs
// fixed step sizes and orders

bool Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), order);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0,i,i);
			invS.set(1,i,i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode, rangeDim+1, i);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, order - 1);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);

		intDifferences.push_back(intTemp);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return false;
	}
	else
	{
		for(int i=0; i<rangeDim; ++i)
			x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		// add the uncertainties and the cutoff intervals onto the result
		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}

bool Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=globalMaxOrder; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode, rangeDim+1, i);	// compute c(t)
	}

	c.nctrunc(orders);

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, orders[i] - 1);
		Adrdt.tms.push_back(tmTemp);
	}
	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return false;
	}
	else
	{
		for(int i=0; i<rangeDim; ++i)
			x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}

// adaptive step sizes and fixed orders

bool Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), order);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode, rangeDim+1, i);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, order - 1);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
	}

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return false;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*order);

		for(int i=0; i<rangeDim; ++i)	// update the uncertainties in the step
		{
			step_uncertainties[i] = step_exp_table[1] * uncertainties[i];
		}

		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}

bool Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=globalMaxOrder; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode, rangeDim+1, i);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, orders[i] - 1);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);
	}

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return false;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*globalMaxOrder);

		for(int i=0; i<rangeDim; ++i)	// update the uncertainties in the step
		{
			step_uncertainties[i] = step_exp_table[1] * uncertainties[i];
		}

		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders);

		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}

// adaptive orders and fixed step sizes

bool Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), order);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode, rangeDim+1, i);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, order - 1);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	// add the uncertainties onto the reault
	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	int newOrder = order;

	for(; !bfound;)
	{
		++newOrder;

		if(newOrder > maxOrder)
		{
			return false;
		}

		// increase the approximation orders by 1
		x.Picard_no_remainder_assign(x0, ode, rangeDim+1, newOrder);

		for(int i=0; i<rangeDim; ++i)	// apply the estimation again
		{
			x.tms[i].remainder = estimation[i] + step_uncertainties[i];
		}

		// compute the Picard operation again
		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, newOrder);

		// Update the irreducible part
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial polyTemp;
			polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

			Interval intTemp;
			polyTemp.intEvalNormal(intTemp, step_exp_table);
			intDifferences[i] = intTemp;
		}

		// add the uncertainties onto the result
		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;
		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		// add the uncertainties
		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	order = newOrder;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}

bool Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const int localMaxOrder, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=localMaxOrder; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode, rangeDim+1, i);	// compute c(t)
	}

	c.nctrunc(orders);

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, orders[i] - 1);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	vector<bool> bIncrease;
	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bIncrease.push_back(true);
			if(bfound)
				bfound = false;
		}
		else
		{
			bIncrease.push_back(false);
		}
	}

	vector<int> newOrders = orders;
	bool bIncreaseOthers = false;
	int numIncrease = 0;

	vector<bool> bIncreased;
	for(int i=0; i<rangeDim; ++i)
	{
		bIncreased.push_back(false);
	}

	for(; !bfound;)
	{
		bool bChanged = false;

		if(bIncreaseOthers)
		{
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(!bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			bIncreaseOthers = false;
		}
		else
		{
			numIncrease = 0;
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					++numIncrease;

					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			if(numIncrease < newOrders.size())
				bIncreaseOthers = true;
		}

		if(!bChanged)
		{
			return false;
		}

		// increase the approximation orders
		x.Picard_no_remainder_assign(x0, ode, rangeDim+1, newOrders, bIncreased);

		for(int i=0; i<rangeDim; ++i)	// apply the estimation again
		{
			x.tms[i].remainder = estimation[i] + step_uncertainties[i];
		}

		// compute the Picard operation again
		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, newOrders);

		for(int i=0; i<rangeDim; ++i)
		{
			// update the irreducible part if necessary
			if(bIncreased[i])
			{
				Polynomial polyTemp;
				polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

				Interval intTemp;
				polyTemp.intEvalNormal(intTemp, step_exp_table);
				intDifferences[i] = intTemp;
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = true;
					if(bfound)
						bfound = false;
				}
				else
				{
					bfound = false;
					break;
				}
			}
			else
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = false;
				}
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			bIncreased[i] = false;
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	orders = newOrders;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return true;
}















// integration scheme for non-polynomial ODEs (using Taylor approximations)
// fixed step sizes and orders

bool Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const
{
	int rangeDim = strOde.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);

	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), order);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde, i, uncertainty_centers);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.order = order - 1;
	parseSetting.flowpipe = c_plus_Ar;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde[i] + suffix;

		parseODE();		// call the parser

		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];	// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, order, uncertainty_centers);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);

		intDifferences.push_back(intTemp);
	}

	// add the uncertainties and the rounded intervals onto the result
	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return false;
	}
	else
	{
		for(int i=0; i<rangeDim; ++i)
			x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], order);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders.push_back(tmvTemp.tms[i].remainder);
		}

		// add the uncertainties and the rounded intervals onto the result
		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	return true;
}

bool Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const
{
	int rangeDim = strOde.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);

	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=globalMaxOrder; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde, i, uncertainty_centers);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.flowpipe = c_plus_Ar;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde[i] + suffix;
		parseSetting.order = orders[i] - 1;

		parseODE();		// call the parser

		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];	// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, orders, uncertainty_centers);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);

		intDifferences.push_back(intTemp);
	}

	// add the uncertainties and the rounded intervals onto the result
	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return false;
	}
	else
	{
		for(int i=0; i<rangeDim; ++i)
			x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], orders);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders.push_back(tmvTemp.tms[i].remainder);
		}

		// add the uncertainties and the rounded intervals onto the result
		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	return true;
}


// adaptive step sizes and fixed orders
bool Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const
{
	int rangeDim = strOde.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);

	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), order);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde, i, uncertainty_centers);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.flowpipe = c_plus_Ar;
	parseSetting.order = order - 1;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde[i] + suffix;

		parseODE();		// call the parser

		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
	}

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, order, uncertainty_centers);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return false;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*order);

		for(int i=0; i<rangeDim; ++i)	// update the uncertainties in the step
		{
			step_uncertainties[i] = step_exp_table[1] * uncertainties[i];
		}

		x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, order, uncertainty_centers);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], order);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	return true;
}

bool Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const
{
	int rangeDim = strOde.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);

	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=globalMaxOrder; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde, i, uncertainty_centers);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.flowpipe = c_plus_Ar;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde[i] + suffix;
		parseSetting.order = orders[i] - 1;

		parseODE();		// call the parser

		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);
	}

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, orders, uncertainty_centers);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return false;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*globalMaxOrder);

		for(int i=0; i<rangeDim; ++i)	// update the uncertainties in the step
		{
			step_uncertainties[i] = step_exp_table[1] * uncertainties[i];
		}

		x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, orders, uncertainty_centers);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], orders);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	return true;
}


// adaptive orders and fixed step sizes
bool Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const
{
	int rangeDim = strOde.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);

	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), order);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde, i, uncertainty_centers);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.order = order - 1;
	parseSetting.flowpipe = c_plus_Ar;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde[i] + suffix;

		parseODE();		// call the parser

		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];	// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, order, uncertainty_centers);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);

		intDifferences.push_back(intTemp);
	}

	// add the uncertainties and the rounded intervals onto the result
	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	int newOrder = order;

	for(; !bfound;)
	{
		++newOrder;

		if(newOrder > maxOrder)
		{
			return false;
		}

		// increase the approximation orders by 1
		x.Picard_non_polynomial_taylor_no_remainder_assign(x0, strOde, newOrder, uncertainty_centers);

		for(int i=0; i<rangeDim; ++i)	// apply the estimation again
		{
			x.tms[i].remainder = estimation[i] + step_uncertainties[i];
		}

		// compute the Picard operation again
		x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, newOrder, uncertainty_centers);

		// Update the irreducible part
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial polyTemp;
			polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

			Interval intTemp;
			polyTemp.intEvalNormal(intTemp, step_exp_table);
			intDifferences[i] = intTemp;
		}

		// add the uncertainties onto the result
		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;
		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], newOrder);

		// add the uncertainties
		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	order = newOrder;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	return true;
}

bool Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const int localMaxOrder, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const
{
	int rangeDim = strOde.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);

	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;
		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);
		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders);

	// we normalized the interval bound of r0
	vector<Interval> boundOfr0;
	result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0, i, i);
			invS.set(1, i, i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table);

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=localMaxOrder; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde, i, uncertainty_centers);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.flowpipe = c_plus_Ar;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde[i] + suffix;
		parseSetting.order = orders[i] - 1;

		parseODE();		// call the parser

		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;
	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders);
	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff();

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	vector<Interval> step_uncertainties;
	for(int i=0; i<rangeDim; ++i)
	{
		step_uncertainties.push_back(step_exp_table[1] * uncertainties[i]);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i] + step_uncertainties[i];	// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, orders, uncertainty_centers);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);

		intDifferences.push_back(intTemp);
	}

	// add the uncertainties and the rounded intervals onto the result
	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += step_uncertainties[i];
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	vector<bool> bIncrease;
	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bIncrease.push_back(true);
			if(bfound)
				bfound = false;
		}
		else
		{
			bIncrease.push_back(false);
		}
	}

	vector<int> newOrders = orders;
	bool bIncreaseOthers = false;
	int numIncrease = 0;

	vector<bool> bIncreased;
	for(int i=0; i<rangeDim; ++i)
	{
		bIncreased.push_back(false);
	}

	for(; !bfound;)
	{
		bool bChanged = false;

		if(bIncreaseOthers)
		{
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(!bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			bIncreaseOthers = false;
		}
		else
		{
			numIncrease = 0;
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					++numIncrease;

					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			if(numIncrease < newOrders.size())
				bIncreaseOthers = true;
		}

		if(!bChanged)
		{
			return false;
		}

		// increase the approximation orders
		x.Picard_non_polynomial_taylor_no_remainder_assign(x0, strOde, newOrders, bIncreased, uncertainty_centers);

		for(int i=0; i<rangeDim; ++i)	// apply the estimation again
		{
			x.tms[i].remainder = estimation[i] + step_uncertainties[i];
		}

		// compute the Picard operation again
		x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, newOrders, uncertainty_centers);

		for(int i=0; i<rangeDim; ++i)
		{
			// update the irreducible part if necessary
			if(bIncreased[i])
			{
				Polynomial polyTemp;
				polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

				Interval intTemp;
				polyTemp.intEvalNormal(intTemp, step_exp_table);
				intDifferences[i] = intTemp;
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += step_uncertainties[i];
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = true;
					if(bfound)
						bfound = false;
				}
				else
				{
					bfound = false;
					break;
				}
			}
			else
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = false;
				}
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			bIncreased[i] = false;
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(; !bfinished;)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], newOrders);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += step_uncertainties[i];
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{
			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	orders = newOrders;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	return true;
}

Flowpipe & Flowpipe::operator = (const Flowpipe & flowpipe)
{
	if(this == &flowpipe)
		return *this;

	tmvPre = flowpipe.tmvPre;
	tmv = flowpipe.tmv;
	domain = flowpipe.domain;
	return *this;
}








































// class Continuous_system

ContinuousSystem::ContinuousSystem()
{
}

ContinuousSystem::ContinuousSystem(const TaylorModelVec & ode_input, const vector<Interval> & uncertainties_input, const Flowpipe & initialSet_input)
{
	int rangeDim = ode_input.tms.size();
	Interval intZero;

	initialSet = initialSet_input;
	tmvOde = ode_input;
	uncertainties = uncertainties_input;

	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		uncertainties[i].remove_midpoint(M);

		if(!M.subseteq(intZero))
		{
			TaylorModel tmTemp(M, rangeDim+1);
			tmvOde.tms[i].add_assign(tmTemp);
		}

		HornerForm hf;
		tmvOde.tms[i].expansion.toHornerForm(hf);
		hfOde.push_back(hf);
	}
}

ContinuousSystem::ContinuousSystem(const vector<string> & strOde_input, const vector<Interval> & uncertainties_input, const Flowpipe & initialSet_input)
{
	uncertainties = uncertainties_input;

	for(int i=0; i<uncertainties.size(); ++i)
	{
		Interval M;
		uncertainties[i].remove_midpoint(M);
		uncertainty_centers.push_back(M);
	}

	strOde = strOde_input;
	initialSet = initialSet_input;
}

ContinuousSystem::ContinuousSystem(const ContinuousSystem & system)
{
	tmvOde				=	system.tmvOde;
	hfOde				=	system.hfOde;
	initialSet			=	system.initialSet;
	uncertainties		=	system.uncertainties;
	uncertainty_centers	=	system.uncertainty_centers;
	strOde				=	system.strOde;
}

ContinuousSystem::~ContinuousSystem()
{
	hfOde.clear();
	uncertainties.clear();
	uncertainty_centers.clear();
	strOde.clear();
}

// fixed step sizes and orders

void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, order);

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, order, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", order);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, orders);

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, orders, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive step sizes and fixed orders

void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, order);

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, order, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("order = %d\n", order);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;

			double tDiffer = time - t;

			if(newStep > tDiffer)
			{
				newStep = tDiffer;
			}

			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, orders);

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, orders, globalMaxOrder, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive orders and fixed step sizes

void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	int newOrder = order;
	int currentMaxOrder = order;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansionHF;
	vector<Polynomial> taylorExpansionMF;
	vector<Polynomial> highestTerms;

	computeTaylorExpansion(taylorExpansionHF, taylorExpansionMF, highestTerms, polyODE, order);

	vector<vector<HornerForm> > expansions;
	expansions.push_back(taylorExpansionHF);

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, expansions[newOrder-order], precondition, step_exp_table, step_end_exp_table, newOrder, maxOrder, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", newOrder);
			}

			if(newOrder > order)
			{
				--newOrder;

				if(newOrder > currentMaxOrder)
				{
					for(int i=currentMaxOrder; i<newOrder; ++i)
					{
						vector<HornerForm> newTaylorExpansionHF;
						vector<Polynomial> newTaylorExpansionMF;

						increaseExpansionOrder(newTaylorExpansionHF, newTaylorExpansionMF, highestTerms, taylorExpansionMF, polyODE, i);

						expansions.push_back(newTaylorExpansionHF);
						taylorExpansionMF = newTaylorExpansionMF;
					}

					currentMaxOrder = newOrder;
				}
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<int> newOrders = orders;
	vector<int> localMaxOrders = orders;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansionHF;
	vector<Polynomial> taylorExpansionMF;
	vector<Polynomial> highestTerms;

	computeTaylorExpansion(taylorExpansionHF, taylorExpansionMF, highestTerms, polyODE, orders);

	vector<vector<HornerForm> > expansions;
	vector<HornerForm> emptySet;
	for(int i=0; i<taylorExpansionHF.size(); ++i)
	{
		expansions.push_back(emptySet);
		expansions[i].push_back(taylorExpansionHF[i]);
	}

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, taylorExpansionHF, precondition, step_exp_table, step_end_exp_table, newOrders, maxOrders, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = newOrders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), newOrders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), newOrders[num]);
			}

			for(int i=0; i<newOrders.size(); ++i)
			{
				if(newOrders[i] > orders[i])
				{
					--newOrders[i];

					if(newOrders[i] > localMaxOrders[i])
					{
						for(int j=localMaxOrders[i]; j<newOrders[i]; ++j)
						{
							HornerForm newTaylorExpansionHF;
							Polynomial newTaylorExpansionMF;

							increaseExpansionOrder(newTaylorExpansionHF, newTaylorExpansionMF, highestTerms[i], taylorExpansionMF[i], polyODE, j);

							expansions[i].push_back(newTaylorExpansionHF);
							taylorExpansionMF[i] = newTaylorExpansionMF;
						}
					}

					localMaxOrders[i] = newOrders[i];

					taylorExpansionHF[i] = expansions[i][newOrders[i]-orders[i]];
				}
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}











// for high-degree ODEs
// fixed step sizes and orders

void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, precondition, step_exp_table, step_end_exp_table, order, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", order);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, precondition, step_exp_table, step_end_exp_table, orders, globalMaxOrder, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive step sizes and fixed orders

void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, order, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("order = %d\n", order);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, orders, globalMaxOrder, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive orders and fixed step sizes

void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	int newOrder = order;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, precondition, step_exp_table, step_end_exp_table, newOrder, maxOrder, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", newOrder);
			}

			if(newOrder > order)
			{
				--newOrder;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	vector<int> newOrders = orders;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int localMaxOrder = newOrders[0];
		for(int i=1; i<newOrders.size(); ++i)
		{
			if(localMaxOrder < newOrders[i])
				localMaxOrder = newOrders[i];
		}

		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, precondition, step_exp_table, step_end_exp_table, newOrders, localMaxOrder, maxOrders, estimation, uncertainties);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = newOrders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), newOrders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), newOrders[num]);
			}

			for(int i=0; i<newOrders.size(); ++i)
			{
				if(newOrders[i] > orders[i])
					--newOrders[i];
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}













// for non-polynomial ODEs (using Taylor approximations)
// fixed step sizes and orders

void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, precondition, step_exp_table, step_end_exp_table, order, estimation, uncertainties, uncertainty_centers);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", order);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, precondition, step_exp_table, step_end_exp_table, orders, globalMaxOrder, estimation, uncertainties, uncertainty_centers);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive step sizes and fixed orders
void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, order, estimation, uncertainties, uncertainty_centers);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("order = %d\n", order);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, orders, globalMaxOrder, estimation, uncertainties, uncertainty_centers);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}


// adaptive orders and fixed step sizes
void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	int newOrder = order;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, precondition, step_exp_table, step_end_exp_table, newOrder, maxOrder, estimation, uncertainties, uncertainty_centers);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", newOrder);
			}

			if(newOrder > order)
			{
				--newOrder;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	vector<int> newOrders = orders;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int localMaxOrder = newOrders[0];
		for(int i=1; i<newOrders.size(); ++i)
		{
			if(localMaxOrder < newOrders[i])
				localMaxOrder = newOrders[i];
		}

		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, precondition, step_exp_table, step_end_exp_table, newOrders, localMaxOrder, maxOrders, estimation, uncertainties, uncertainty_centers);

		if(bvalid)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = newOrders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), newOrders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), newOrders[num]);
			}

			for(int i=0; i<newOrders.size(); ++i)
			{
				if(newOrders[i] > orders[i])
					--newOrders[i];
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

ContinuousSystem & ContinuousSystem::operator = (const ContinuousSystem & system)
{
	if(this == &system)
		return *this;

	tmvOde				=	system.tmvOde;
	hfOde				=	system.hfOde;
	initialSet			=	system.initialSet;
	uncertainties		=	system.uncertainties;
	uncertainty_centers	=	system.uncertainty_centers;
	strOde				=	system.strOde;

	return *this;
}




































// class ContinuousReachability

ContinuousReachability::ContinuousReachability()
{
}

ContinuousReachability::~ContinuousReachability()
{
	outputAxes.clear();
	flowpipes.clear();
	orders.clear();
	maxOrders.clear();
	flowpipesCompo.clear();
	domains.clear();
	unsafeSet.clear();
	stateVarTab.clear();
	stateVarNames.clear();
	tmVarTab.clear();
	tmVarNames.clear();
}

void ContinuousReachability::dump(FILE *fp) const
{
	fprintf(fp,"state var ");
	for(int i=0; i<stateVarNames.size()-1; ++i)
	{
		fprintf(fp, "%s,", stateVarNames[i].c_str());
	}
	fprintf(fp, "%s\n\n", stateVarNames[stateVarNames.size()-1].c_str());

	switch(plotFormat)
	{
	case PLOT_GNUPLOT:
		switch(plotSetting)
		{
		case PLOT_INTERVAL:
			fprintf(fp, "gnuplot interval %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_OCTAGON:
			fprintf(fp, "gnuplot octagon %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_GRID:
			fprintf(fp, "gnuplot grid %d %s , %s\n\n", numSections, stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		}
		break;
	case PLOT_MATLAB:
		switch(plotSetting)
		{
		case PLOT_INTERVAL:
			fprintf(fp, "matlab interval %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_OCTAGON:
			fprintf(fp, "matlab octagon %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_GRID:
			fprintf(fp, "matlab grid %d %s , %s\n\n", numSections, stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		}
		break;
	}

	fprintf(fp, "output %s\n\n", outputFileName);

	if(bSafetyChecking)
	{
		// dump the unsafe set
		fprintf(fp, "unsafe set\n{\n");

		for(int i=0; i<unsafeSet.size(); ++i)
		{
			unsafeSet[i].dump(fp, stateVarNames);
		}

		fprintf(fp, "}\n\n");
	}

	fprintf(fp, "continuous flowpipes\n{\n");

	fprintf(fp, "tm var ");
	for(int i=1; i<tmVarNames.size()-1; ++i)
	{
		fprintf(fp, "%s,", tmVarNames[i].c_str());
	}
	fprintf(fp, "%s\n\n", tmVarNames[tmVarNames.size()-1].c_str());

	list<TaylorModelVec>::const_iterator fpIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	// every Taylor model flowpipe is enclosed by braces

	for(; fpIter != flowpipesCompo.end(); ++fpIter, ++doIter)
	{
		fprintf(fp, "{\n");
		fpIter->dump_interval(fp, stateVarNames, tmVarNames);

		//dump the domain
		for(int i=0; i<doIter->size(); ++i)
		{
			fprintf(fp, "%s in ", tmVarNames[i].c_str());
			(*doIter)[i].dump(fp);
			fprintf(fp, "\n");
		}

		fprintf(fp, "}\n\n");
	}

	fprintf(fp, "}\n");
}

void ContinuousReachability::run()
{
	compute_factorial_rec(globalMaxOrder+1);
	compute_power_4(globalMaxOrder+1);
	compute_double_factorial(2*globalMaxOrder);

	switch(integrationScheme)
	{
	case LOW_DEGREE:
	{
		switch(orderType)
		{
		case UNIFORM:
			if(bAdaptiveSteps)
			{
				system.reach_low_degree(flowpipes, step, miniStep, time, orders[0], precondition, estimation, bPrint, stateVarNames);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_low_degree(flowpipes, step, time, orders[0], maxOrders[0], precondition, estimation, bPrint, stateVarNames);
			}
			else
			{
				system.reach_low_degree(flowpipes, step, time, orders[0], precondition, estimation, bPrint, stateVarNames);
			}
			break;
		case MULTI:
			if(bAdaptiveSteps)
			{
				system.reach_low_degree(flowpipes, step, miniStep, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_low_degree(flowpipes, step, time, orders, maxOrders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames);
			}
			else
			{
				system.reach_low_degree(flowpipes, step, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames);
			}
			break;
		}
		break;
	}

	case HIGH_DEGREE:
	{
		switch(orderType)
		{
		case UNIFORM:
			if(bAdaptiveSteps)
			{
				system.reach_high_degree(flowpipes, step, miniStep, time, orders[0], precondition, estimation, bPrint, stateVarNames);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_high_degree(flowpipes, step, time, orders[0], maxOrders[0], precondition, estimation, bPrint, stateVarNames);
			}
			else
			{
				system.reach_high_degree(flowpipes, step, time, orders[0], precondition, estimation, bPrint, stateVarNames);
			}
			break;
		case MULTI:
			if(bAdaptiveSteps)
			{
				system.reach_high_degree(flowpipes, step, miniStep, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_high_degree(flowpipes, step, time, orders, maxOrders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames);
			}
			else
			{
				system.reach_high_degree(flowpipes, step, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames);
			}
			break;
		}
		break;
	}

	case NONPOLY_TAYLOR:
	{
		switch(orderType)
		{
		case UNIFORM:
			if(bAdaptiveSteps)
			{
				system.reach_non_polynomial_taylor(flowpipes, step, miniStep, time, orders[0], precondition, estimation, bPrint, stateVarNames);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_non_polynomial_taylor(flowpipes, step, time, orders[0], maxOrders[0], precondition, estimation, bPrint, stateVarNames);
			}
			else
			{
				system.reach_non_polynomial_taylor(flowpipes, step, time, orders[0], precondition, estimation, bPrint, stateVarNames);
			}
			break;
		case MULTI:
			if(bAdaptiveSteps)
			{
				system.reach_non_polynomial_taylor(flowpipes, step, miniStep, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_non_polynomial_taylor(flowpipes, step, time, orders, maxOrders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames);
			}
			else
			{
				system.reach_non_polynomial_taylor(flowpipes, step, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames);
			}
			break;
		}
		break;
	}
	}
}

void ContinuousReachability::composition()
{
	flowpipesCompo.clear();
	domains.clear();

	vector<Interval> step_exp_table;
	Interval intStep;

	list<Flowpipe>::const_iterator iter;

	for(iter = flowpipes.begin(); iter != flowpipes.end(); ++iter)
	{
		if(step_exp_table.size() == 0 || intStep != iter->domain[0])
		{
			construct_step_exp_table(step_exp_table, iter->domain[0], globalMaxOrder);
			intStep = iter->domain[0];
		}

		TaylorModelVec tmvTemp;

		iter->composition_normal(tmvTemp, step_exp_table);

		flowpipesCompo.push_back(tmvTemp);
		domains.push_back(iter->domain);
	}
}

int ContinuousReachability::safetyChecking() const
{
	if(unsafeSet.size() == 0)
	{
		return UNSAFE;	// since the whole state space is unsafe, the system is not safe
	}

	bool bDumpCounterexamples = true;

	int mkres = mkdir(counterexampleDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for counterexamples.\n");
		bDumpCounterexamples = false;
	}

	char filename_counterexamples[NAME_SIZE+10];
	FILE *fpDumpCounterexamples;

	if(bDumpCounterexamples)
	{
		sprintf(filename_counterexamples, "%s%s%s", counterexampleDir, outputFileName, str_counterexample_dumping_name_suffix);
		fpDumpCounterexamples = fopen(filename_counterexamples, "w");
	}

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	vector<Interval> step_exp_table;

	int rangeDim = tmvIter->tms.size();
	int domainDim = doIter->size();

	int result = SAFE;

	list<TaylorModelVec> unsafe_flowpipes;
	list<vector<Interval> > unsafe_flowpipe_domains;
	list<Interval> globalTimes;

	int maxOrder = 0;
	Interval globalTime;

	for(; tmvIter!=flowpipesCompo.end(); ++tmvIter, ++doIter)
	{
		int tmp = maxOrder;
		for(int i=0; i<tmvIter->tms.size(); ++i)
		{
			int order = tmvIter->tms[i].expansion.degree();
			if(maxOrder < order)
			{
				maxOrder = order;
			}
		}

		if(step_exp_table.size() == 0 || step_exp_table[1] != (*doIter)[0] || maxOrder > tmp)
		{
			construct_step_exp_table(step_exp_table, (*doIter)[0], 2*maxOrder);
		}

		bool bsafe = false;

		vector<Interval> tmvPolyRange;
		tmvIter->polyRangeNormal(tmvPolyRange, step_exp_table);

		for(int i=0; i<unsafeSet.size(); ++i)
		{
			TaylorModel tmTemp;

			// interval evaluation on the constraint
			unsafeSet[i].hf.insert_normal(tmTemp, *tmvIter, tmvPolyRange, step_exp_table, domainDim);

			Interval intTemp;
			tmTemp.intEvalNormal(intTemp, step_exp_table);

			if(intTemp > unsafeSet[i].B)
			{
				// no intersection with the unsafe set
				bsafe = true;
				break;
			}
			else
			{
				continue;
			}
		}

		if(!bsafe)
		{
			unsafe_flowpipes.push_back(*tmvIter);
			unsafe_flowpipe_domains.push_back(*doIter);
			globalTimes.push_back(globalTime);

			result = UNKNOWN;
		}

		globalTime += (*doIter)[0];
	}

	if(bDumpCounterexamples && unsafe_flowpipes.size() > 0)
	{
		dump_potential_counterexample(fpDumpCounterexamples, unsafe_flowpipes, unsafe_flowpipe_domains, globalTimes);
	}

	if(bDumpCounterexamples)
	{
		fclose(fpDumpCounterexamples);
	}

	return result;
}

unsigned long ContinuousReachability::numOfFlowpipes() const
{
	return (unsigned long)(flowpipes.size()) - 1;
}

void ContinuousReachability::dump_potential_counterexample(FILE *fp, const list<TaylorModelVec> & flowpipes, const list<vector<Interval> > & domains, const list<Interval> & globalTimes) const
{
	list<TaylorModelVec>::const_iterator fpIter = flowpipes.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();
	list<Interval>::const_iterator timeIter = globalTimes.begin();

	for(; fpIter!=flowpipes.end(); ++fpIter, ++doIter, ++timeIter)
	{
		fprintf(fp, "starting time %lf\n{\n", timeIter->sup());

		fpIter->dump_interval(fp, stateVarNames, tmVarNames);

		for(int i=0; i<doIter->size(); ++i)
		{
			fprintf(fp, "%s in ", tmVarNames[i].c_str());
			(*doIter)[i].dump(fp);
			fprintf(fp, "\n");
		}

		fprintf(fp, "}\n\n\n");
	}
}

void ContinuousReachability::plot_2D() const
{
	char filename[NAME_SIZE+10];

	switch(plotFormat)
	{
	case PLOT_GNUPLOT:
		sprintf(filename, "%s%s.plt", outputDir, outputFileName);
		break;
	case PLOT_MATLAB:
		sprintf(filename, "%s%s.m", outputDir, outputFileName);
		break;
	}

	FILE *fpPlotting = fopen(filename, "w");

	if(fpPlotting == NULL)
	{
		printf("Can not create the plotting file.\n");
		exit(1);
	}

	printf("Generating the plotting file...\n");
	switch(plotFormat)
	{
	case PLOT_GNUPLOT:
		plot_2D_GNUPLOT(fpPlotting);
		break;
	case PLOT_MATLAB:
		plot_2D_MATLAB(fpPlotting);
		break;
	}
	printf("Done.\n");

	fclose(fpPlotting);
}

void ContinuousReachability::plot_2D_GNUPLOT(FILE *fp) const
{
	switch(plotSetting)
	{
	case PLOT_INTERVAL:
		plot_2D_interval_GNUPLOT(fp);
		break;
	case PLOT_OCTAGON:
		plot_2D_octagon_GNUPLOT(fp);
		break;
	case PLOT_GRID:
		plot_2D_grid_GNUPLOT(fp);
		break;
	}
}

void ContinuousReachability::plot_2D_interval_GNUPLOT(FILE *fp) const
{
	fprintf(fp, "set terminal postscript\n");

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.eps", imageDir, outputFileName);
	fprintf(fp, "set output '%s'\n", filename);

	fprintf(fp, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(fp, "set autoscale\n");
	fprintf(fp, "unset label\n");
	fprintf(fp, "set xtic auto\n");
	fprintf(fp, "set ytic auto\n");
	fprintf(fp, "set xlabel \"%s\"\n", stateVarNames[outputAxes[0]].c_str());
	fprintf(fp, "set ylabel \"%s\"\n", stateVarNames[outputAxes[1]].c_str());
	fprintf(fp, "plot '-' notitle with lines ls 1\n");

	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intStep;

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	int maxOrder = 0;

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		int tmp = maxOrder;
		for(int i=0; i<tmvIter->tms.size(); ++i)
		{
			int order = tmvIter->tms[i].expansion.degree();
			if(maxOrder < order)
			{
				maxOrder = order;
			}
		}

		if(step_exp_table.size() == 0 || intStep != (*doIter)[0] || maxOrder > tmp)
		{
			construct_step_exp_table(step_exp_table, step_end_exp_table, (*doIter)[0].sup(), maxOrder);
			intStep = (*doIter)[0];
		}

		vector<Interval> box;
		tmvIter->intEvalNormal(box, step_exp_table);

		Interval X(box[outputAxes[0]]), Y(box[outputAxes[1]]);

		// output the vertices
		fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
		fprintf(fp, "%lf %lf\n", X.sup(), Y.inf());
		fprintf(fp, "%lf %lf\n", X.sup(), Y.sup());
		fprintf(fp, "%lf %lf\n", X.inf(), Y.sup());
		fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
		fprintf(fp, "\n\n");
	}

	fprintf(fp, "e\n");
}

void ContinuousReachability::plot_2D_octagon_GNUPLOT(FILE *fp) const
{
	int x = outputAxes[0];
	int y = outputAxes[1];

	int rangeDim = stateVarNames.size();
	Matrix output_poly_temp(8, rangeDim);

	output_poly_temp.set(1, 0, x);
	output_poly_temp.set(1, 1, y);
	output_poly_temp.set(-1, 2, x);
	output_poly_temp.set(-1, 3, y);
	output_poly_temp.set(1/sqrt(2), 4, x);
	output_poly_temp.set(1/sqrt(2), 4, y);
	output_poly_temp.set(1/sqrt(2), 5, x);
	output_poly_temp.set(-1/sqrt(2), 5, y);
	output_poly_temp.set(-1/sqrt(2), 6, x);
	output_poly_temp.set(1/sqrt(2), 6, y);
	output_poly_temp.set(-1/sqrt(2), 7, x);
	output_poly_temp.set(-1/sqrt(2), 7, y);

	// Construct the 2D template matrix.
	int rows = 8;
	int cols = rangeDim;

	Matrix sortedTemplate(rows, cols);
	RowVector rowVec(cols);
	list<RowVector> sortedRows;
	list<RowVector>::iterator iterp, iterq;

	output_poly_temp.row(rowVec, 0);
	sortedRows.push_back(rowVec);

	bool bInserted;

	// Sort the row vectors in the template by anti-clockwise order (only in the x-y space).
	for(int i=1; i<rows; ++i)
	{
		iterp = sortedRows.begin();
		iterq = iterp;
		++iterq;
		bInserted = false;

		for(; iterq != sortedRows.end();)
		{
			double tmp1 = output_poly_temp.get(i,x) * iterp->get(y) - output_poly_temp.get(i,y) * iterp->get(x);
			double tmp2 = output_poly_temp.get(i,x) * iterq->get(y) - output_poly_temp.get(i,y) * iterq->get(x);

			if(tmp1 < 0 && tmp2 > 0)
			{
				output_poly_temp.row(rowVec, i);
				sortedRows.insert(iterq, rowVec);
				bInserted = true;
				break;
			}
			else
			{
				++iterp;
				++iterq;
			}
		}

		if(!bInserted)
		{
			output_poly_temp.row(rowVec, i);
			sortedRows.push_back(rowVec);
		}
	}

	iterp = sortedRows.begin();
	for(int i=0; i<rows; ++i, ++iterp)
	{
		for(int j=0; j<cols; ++j)
		{
			sortedTemplate.set(iterp->get(j), i, j);
		}
	}

	ColVector b(rows);
	Polyhedron polyTemplate(sortedTemplate, b);

	fprintf(fp, "set terminal postscript\n");

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.eps", imageDir, outputFileName);
	fprintf(fp, "set output '%s'\n", filename);

	fprintf(fp, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(fp, "set autoscale\n");
	fprintf(fp, "unset label\n");
	fprintf(fp, "set xtic auto\n");
	fprintf(fp, "set ytic auto\n");
	fprintf(fp, "set xlabel \"%s\"\n", stateVarNames[outputAxes[0]].c_str());
	fprintf(fp, "set ylabel \"%s\"\n", stateVarNames[outputAxes[1]].c_str());
	fprintf(fp, "plot '-' notitle with lines ls 1\n");

	// Compute the intersections of two facets.
	// The vertices are ordered clockwisely.

	gsl_matrix *C = gsl_matrix_alloc(2,2);
	gsl_vector *d = gsl_vector_alloc(2);
	gsl_vector *vertex = gsl_vector_alloc(2);

	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intStep;

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	int maxOrder = 0;

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		int tmp = maxOrder;
		for(int i=0; i<tmvIter->tms.size(); ++i)
		{
			int order = tmvIter->tms[i].expansion.degree();
			if(maxOrder < order)
			{
				maxOrder = order;
			}
		}

		if(step_exp_table.size() == 0 || intStep != (*doIter)[0] || maxOrder > tmp)
		{
			construct_step_exp_table(step_exp_table, step_end_exp_table, (*doIter)[0].sup(), maxOrder);
			intStep = (*doIter)[0];
		}


		templatePolyhedronNormal(polyTemplate, *tmvIter, step_exp_table);

		double f1, f2;

		list<LinearConstraint>::iterator iterp, iterq;
		iterp = iterq = polyTemplate.constraints.begin();
		++iterq;

		for(; iterq != polyTemplate.constraints.end(); ++iterp, ++iterq)
		{
			gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
			gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
			gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
			gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

			gsl_vector_set(d, 0, iterp->B.midpoint());
			gsl_vector_set(d, 1, iterq->B.midpoint());

			gsl_linalg_HH_solve(C, d, vertex);

			double v1 = gsl_vector_get(vertex, 0);
			double v2 = gsl_vector_get(vertex, 1);

			if(iterp == polyTemplate.constraints.begin())
			{
				f1 = v1;
				f2 = v2;
			}

			fprintf(fp, "%lf %lf\n", v1, v2);
		}

		iterp = polyTemplate.constraints.begin();
		--iterq;

		gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
		gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
		gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
		gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

		gsl_vector_set(d, 0, iterp->B.midpoint());
		gsl_vector_set(d, 1, iterq->B.midpoint());

		gsl_linalg_HH_solve(C, d, vertex);

		double v1 = gsl_vector_get(vertex, 0);
		double v2 = gsl_vector_get(vertex, 1);

		fprintf(fp, "%lf %lf\n", v1, v2);

		fprintf(fp, "%lf %lf\n", f1, f2);
		fprintf(fp, "\n\n");
	}

	fprintf(fp, "e\n");

	gsl_matrix_free(C);
	gsl_vector_free(d);
	gsl_vector_free(vertex);
}

void ContinuousReachability::plot_2D_grid_GNUPLOT(FILE *fp) const
{
	fprintf(fp, "set terminal postscript\n");

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.eps", imageDir, outputFileName);
	fprintf(fp, "set output '%s'\n", filename);

	fprintf(fp, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(fp, "set autoscale\n");
	fprintf(fp, "unset label\n");
	fprintf(fp, "set xtic auto\n");
	fprintf(fp, "set ytic auto\n");
	fprintf(fp, "set xlabel \"%s\"\n", stateVarNames[outputAxes[0]].c_str());
	fprintf(fp, "set ylabel \"%s\"\n", stateVarNames[outputAxes[1]].c_str());
	fprintf(fp, "plot '-' notitle with lines ls 1\n");

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();
	int domainDim = doIter->size();

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		// decompose the domain
		list<vector<Interval> > grids;

		gridBox(grids, *doIter, numSections);

		// we only consider the output dimensions
		HornerForm hfOutputX;
		Interval remainderX;
		tmvIter->tms[outputAxes[0]].toHornerForm(hfOutputX, remainderX);

		HornerForm hfOutputY;
		Interval remainderY;
		tmvIter->tms[outputAxes[1]].toHornerForm(hfOutputY, remainderY);

		// evaluate the images from all of the grids
		list<vector<Interval> >::const_iterator gIter = grids.begin();
		for(; gIter!=grids.end(); ++gIter)
		{
			Interval X;
			hfOutputX.intEval(X, *gIter);
			X += remainderX;

			Interval Y;
			hfOutputY.intEval(Y, *gIter);
			Y += remainderY;

			// output the vertices
			fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
			fprintf(fp, "%lf %lf\n", X.sup(), Y.inf());
			fprintf(fp, "%lf %lf\n", X.sup(), Y.sup());
			fprintf(fp, "%lf %lf\n", X.inf(), Y.sup());
			fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
			fprintf(fp, "\n\n");
		}
	}

	fprintf(fp, "e\n");
}

void ContinuousReachability::plot_2D_MATLAB(FILE *fp) const
{
	switch(plotSetting)
	{
	case PLOT_INTERVAL:
		plot_2D_interval_MATLAB(fp);
		break;
	case PLOT_OCTAGON:
		plot_2D_octagon_MATLAB(fp);
		break;
	case PLOT_GRID:
		plot_2D_grid_MATLAB(fp);
		break;
	}
}

void ContinuousReachability::plot_2D_interval_MATLAB(FILE *fp) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intStep;

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	int maxOrder = 0;

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		int tmp = maxOrder;
		for(int i=0; i<tmvIter->tms.size(); ++i)
		{
			int order = tmvIter->tms[i].expansion.degree();
			if(maxOrder < order)
			{
				maxOrder = order;
			}
		}

		if(step_exp_table.size() == 0 || intStep != (*doIter)[0] || maxOrder > tmp)
		{
			construct_step_exp_table(step_exp_table, step_end_exp_table, (*doIter)[0].sup(), maxOrder);
			intStep = (*doIter)[0];
		}

		vector<Interval> box;
		tmvIter->intEvalNormal(box, step_exp_table);

		Interval X(box[outputAxes[0]]), Y(box[outputAxes[1]]);

		// output all the vertices
		fprintf(fp,"plot( [%lf,%lf,%lf,%lf,%lf] , [%lf,%lf,%lf,%lf,%lf] , 'b');\nhold on;\nclear;\n",
				X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
	}
}

void ContinuousReachability::plot_2D_octagon_MATLAB(FILE *fp) const
{
	int x = outputAxes[0];
	int y = outputAxes[1];

	int rangeDim = stateVarNames.size();
	Matrix output_poly_temp(8, rangeDim);

	output_poly_temp.set(1, 0, x);
	output_poly_temp.set(1, 1, y);
	output_poly_temp.set(-1, 2, x);
	output_poly_temp.set(-1, 3, y);
	output_poly_temp.set(1/sqrt(2), 4, x);
	output_poly_temp.set(1/sqrt(2), 4, y);
	output_poly_temp.set(1/sqrt(2), 5, x);
	output_poly_temp.set(-1/sqrt(2), 5, y);
	output_poly_temp.set(-1/sqrt(2), 6, x);
	output_poly_temp.set(1/sqrt(2), 6, y);
	output_poly_temp.set(-1/sqrt(2), 7, x);
	output_poly_temp.set(-1/sqrt(2), 7, y);

	// Construct the 2D template matrix.
	int rows = 8;
	int cols = rangeDim;

	Matrix sortedTemplate(rows, cols);
	RowVector rowVec(cols);
	list<RowVector> sortedRows;
	list<RowVector>::iterator iterp, iterq;

	output_poly_temp.row(rowVec, 0);
	sortedRows.push_back(rowVec);

	bool bInserted;

	// Sort the row vectors in the template by anti-clockwise order (only in the x-y space).
	for(int i=1; i<rows; ++i)
	{
		iterp = sortedRows.begin();
		iterq = iterp;
		++iterq;
		bInserted = false;

		for(; iterq != sortedRows.end();)
		{
			double tmp1 = output_poly_temp.get(i,x) * iterp->get(y) - output_poly_temp.get(i,y) * iterp->get(x);
			double tmp2 = output_poly_temp.get(i,x) * iterq->get(y) - output_poly_temp.get(i,y) * iterq->get(x);

			if(tmp1 < 0 && tmp2 > 0)
			{
				output_poly_temp.row(rowVec, i);
				sortedRows.insert(iterq, rowVec);
				bInserted = true;
				break;
			}
			else
			{
				++iterp;
				++iterq;
			}
		}

		if(!bInserted)
		{
			output_poly_temp.row(rowVec, i);
			sortedRows.push_back(rowVec);
		}
	}

	iterp = sortedRows.begin();
	for(int i=0; i<rows; ++i, ++iterp)
	{
		for(int j=0; j<cols; ++j)
		{
			sortedTemplate.set(iterp->get(j), i, j);
		}
	}

	ColVector b(rows);
	Polyhedron polyTemplate(sortedTemplate, b);

	// Compute the intersections of two facets.
	// The vertices are ordered clockwisely.

	gsl_matrix *C = gsl_matrix_alloc(2,2);
	gsl_vector *d = gsl_vector_alloc(2);
	gsl_vector *vertex = gsl_vector_alloc(2);

	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intStep;

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	int maxOrder = 0;

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		int tmp = maxOrder;
		for(int i=0; i<tmvIter->tms.size(); ++i)
		{
			int order = tmvIter->tms[i].expansion.degree();
			if(maxOrder < order)
			{
				maxOrder = order;
			}
		}

		if(step_exp_table.size() == 0 || intStep != (*doIter)[0] || maxOrder > tmp)
		{
			construct_step_exp_table(step_exp_table, step_end_exp_table, (*doIter)[0].sup(), maxOrder);
			intStep = (*doIter)[0];
		}


		templatePolyhedronNormal(polyTemplate, *tmvIter, step_exp_table);

		double f1, f2;

		list<LinearConstraint>::iterator iterp, iterq;
		iterp = iterq = polyTemplate.constraints.begin();
		++iterq;

		vector<double> vertices_x, vertices_y;

		for(; iterq != polyTemplate.constraints.end(); ++iterp, ++iterq)
		{
			gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
			gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
			gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
			gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

			gsl_vector_set(d, 0, iterp->B.midpoint());
			gsl_vector_set(d, 1, iterq->B.midpoint());

			gsl_linalg_HH_solve(C, d, vertex);

			double v1 = gsl_vector_get(vertex, 0);
			double v2 = gsl_vector_get(vertex, 1);

			if(iterp == polyTemplate.constraints.begin())
			{
				f1 = v1;
				f2 = v2;
			}

			vertices_x.push_back(v1);
			vertices_y.push_back(v2);
		}

		iterp = polyTemplate.constraints.begin();
		--iterq;

		gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
		gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
		gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
		gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

		gsl_vector_set(d, 0, iterp->B.midpoint());
		gsl_vector_set(d, 1, iterq->B.midpoint());

		gsl_linalg_HH_solve(C, d, vertex);

		double v1 = gsl_vector_get(vertex, 0);
		double v2 = gsl_vector_get(vertex, 1);

		vertices_x.push_back(v1);
		vertices_y.push_back(v2);
		vertices_x.push_back(f1);
		vertices_y.push_back(f2);

		fprintf(fp, "plot( ");

		fprintf(fp, "[ ");
		for(int i=0; i<vertices_x.size()-1; ++i)
		{
			fprintf(fp, "%lf , ", vertices_x[i]);
		}
		fprintf(fp, "%lf ] , ", vertices_x.back());

		fprintf(fp, "[ ");
		for(int i=0; i<vertices_y.size()-1; ++i)
		{
			fprintf(fp, "%lf , ", vertices_y[i]);
		}
		fprintf(fp, "%lf ] , ", vertices_y.back());

		fprintf(fp, "'b');\nhold on;\nclear;\n");
	}

	gsl_matrix_free(C);
	gsl_vector_free(d);
	gsl_vector_free(vertex);
}

void ContinuousReachability::plot_2D_grid_MATLAB(FILE *fp) const
{
	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();
	int domainDim = doIter->size();

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		// decompose the domain
		list<vector<Interval> > grids;

		gridBox(grids, *doIter, numSections);

		// we only consider the output dimensions
		HornerForm hfOutputX;
		Interval remainderX;
		tmvIter->tms[outputAxes[0]].toHornerForm(hfOutputX, remainderX);

		HornerForm hfOutputY;
		Interval remainderY;
		tmvIter->tms[outputAxes[1]].toHornerForm(hfOutputY, remainderY);

		// evaluate the images from all of the grids
		list<vector<Interval> >::const_iterator gIter = grids.begin();
		for(; gIter!=grids.end(); ++gIter)
		{
			Interval X;
			hfOutputX.intEval(X, *gIter);
			X += remainderX;

			Interval Y;
			hfOutputY.intEval(Y, *gIter);
			Y += remainderY;

			// output the vertices
			fprintf(fp,"plot( [%lf,%lf,%lf,%lf,%lf] , [%lf,%lf,%lf,%lf,%lf] , 'b');\nhold on;\nclear;\n",
					X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
		}
	}
}

bool ContinuousReachability::declareStateVar(const string & vName)
{
	map<string,int>::const_iterator iter;

	if((iter = stateVarTab.find(vName)) == stateVarTab.end())
	{
		stateVarTab[vName] = stateVarNames.size();
		stateVarNames.push_back(vName);
		return true;
	}
	else
	{
		return false;
	}
}

int ContinuousReachability::getIDForStateVar(const string & vName) const
{
	map<string,int>::const_iterator iter;
	if((iter = stateVarTab.find (vName)) == stateVarTab.end())
	{
		return -1;
	}

	return iter -> second;
}

bool ContinuousReachability::getStateVarName(string & vName, int id) const
{
	if(id >= 0 && id<stateVarNames.size())
	{
		vName = stateVarNames[id];
		return true;
	}
	else
	{
		return false;
	}
}

bool ContinuousReachability::declareTMVar(const string & vName)
{
	map<string,int>::const_iterator iter;

	if((iter = tmVarTab.find(vName)) == tmVarTab.end())
	{
		tmVarTab[vName] = tmVarNames.size();
		tmVarNames.push_back(vName);
		return true;
	}
	else
	{
		return false;
	}
}

int ContinuousReachability::getIDForTMVar(const string & vName) const
{
	map<string,int>::const_iterator iter;
	if((iter = tmVarTab.find(vName)) == tmVarTab.end())
	{
		return -1;
	}

	return iter->second;
}

bool ContinuousReachability::getTMVarName(string & vName, int id) const
{
	if(id>=0 && id<tmVarNames.size())
	{
		vName = tmVarNames[id];
		return true;
	}
	else
	{
		return false;
	}
}

void computeTaylorExpansion(TaylorModelVec & result, const TaylorModelVec & first_order_deriv, const TaylorModelVec & ode, const int order)
{
	vector<Interval> intVecZero;
	Interval intZero, intOne(1,1);

	intVecZero.push_back(intOne);
	intVecZero.push_back(intZero);

	// we compute the Taylor expansion (without the 0-order term)
	TaylorModelVec taylorExpansion;
	first_order_deriv.evaluate_t(taylorExpansion, intVecZero);
//	taylorExpansion.nctrunc(order - 1);
	taylorExpansion.mul_assign(0, 1);

	TaylorModelVec tmvLieDeriv_n = first_order_deriv;

	for(int i=2; i<=order; ++i)
	{
		TaylorModelVec tmvTemp;
		tmvLieDeriv_n.LieDerivative_no_remainder(tmvTemp, ode, order - i);

		TaylorModelVec tmvTemp2;
		tmvTemp.evaluate_t(tmvTemp2, intVecZero);
		tmvTemp2.mul_assign(factorial_rec[i]);
		tmvTemp2.mul_assign(0,i);			// multiplied by t^i
//		tmvTemp2.nctrunc(order);

		taylorExpansion.add_assign(tmvTemp2);

		tmvLieDeriv_n = tmvTemp;
	}

	for(int i=0; i<taylorExpansion.tms.size(); ++i)
	{
		taylorExpansion.tms[i].cutoff();
	}

	result = taylorExpansion;
}

void computeTaylorExpansion(TaylorModelVec & result, const TaylorModelVec & first_order_deriv, const TaylorModelVec & ode, const vector<int> & orders)
{
	int rangeDim = ode.tms.size();
	vector<Interval> intVecZero;
	Interval intZero, intOne(1,1);

	intVecZero.push_back(intOne);
	intVecZero.push_back(intZero);

	// we compute the Taylor expansion (without the 0-order term)
	TaylorModelVec taylorExpansion;
	first_order_deriv.evaluate_t(taylorExpansion, intVecZero);
/*
	for(int i=0; i<rangeDim; ++i)
	{
		taylorExpansion.tms[i].nctrunc(orders[i] - 1);
	}
*/
	taylorExpansion.mul_assign(0, 1);

	TaylorModelVec tmvLieDeriv_n = first_order_deriv;

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=2; j<=orders[i]; ++j)
		{
			TaylorModel tmTemp;
			tmvLieDeriv_n.tms[i].LieDerivative_no_remainder(tmTemp, ode, orders[i] - j);

			TaylorModel tmTemp2;
			tmTemp.evaluate_t(tmTemp2, intVecZero);
			tmTemp2.mul_assign(factorial_rec[j]);
			tmTemp2.mul_assign(0,j);
//			tmTemp2.nctrunc(orders[i]);

			taylorExpansion.tms[i].add_assign(tmTemp2);

			tmvLieDeriv_n.tms[i] = tmTemp;
		}
	}

	for(int i=0; i<taylorExpansion.tms.size(); ++i)
	{
		taylorExpansion.tms[i].cutoff();
	}

	result = taylorExpansion;
}

void construct_step_exp_table(vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const int order)
{
	step_exp_table.clear();
	step_end_exp_table.clear();

	Interval intProd(1), intStep(0,step);

	for(int i=0; i<=order; ++i)
	{
		step_exp_table.push_back(intProd);

		Interval intTend(intProd.sup());
		step_end_exp_table.push_back(intTend);

		intProd *= intStep;
	}
}

void construct_step_exp_table(vector<Interval> & step_exp_table, const Interval & step, const int order)
{
	step_exp_table.clear();

	Interval intProd(1);
	step_exp_table.push_back(intProd);

	for(int i=1; i<=order; ++i)
	{
		intProd *= step;
		step_exp_table.push_back(intProd);
	}
}

void preconditionQR(Matrix & result, const TaylorModelVec & x0, const int rangeDim, const int domainDim)
{
	Interval intZero;
	vector<vector<Interval> > intCoefficients;

	for(int i=0; i<rangeDim; ++i)
	{
		vector<Interval> intVecTemp;
		for(int j=0; j<domainDim; ++j)
		{
			intVecTemp.push_back(intZero);
		}
		intCoefficients.push_back(intVecTemp);
	}

	x0.linearCoefficients(intCoefficients);
	Matrix matCoefficients(rangeDim, rangeDim);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=1; j<=rangeDim; ++j)
		{
			matCoefficients.set(intCoefficients[i][j].midpoint(), i, j-1);
		}
	}

	matCoefficients.sortColumns();
	matCoefficients.QRfactor(result);
}

Interval rho(const TaylorModelVec & tmv, const vector<Interval> & l, const vector<Interval> & domain)
{
	int d = l.size();
	TaylorModel tmObj;

	for(int i=0; i<d; ++i)
	{
		TaylorModel tmTemp;
		tmv.tms[i].mul(tmTemp, l[i]);
		tmObj.add_assign(tmTemp);
	}

	Interval intRange;
	tmObj.intEval(intRange, domain);

	Interval S;
	intRange.sup(S);

	return S;
}

Interval rhoNormal(const TaylorModelVec & tmv, const vector<Interval> & l, const vector<Interval> & step_exp_table)
{
	int d = l.size();
	TaylorModel tmObj;

	for(int i=0; i<d; ++i)
	{
		TaylorModel tmTemp;
		tmv.tms[i].mul(tmTemp, l[i]);
		tmObj.add_assign(tmTemp);
	}

	Interval intRange;
	tmObj.intEvalNormal(intRange, step_exp_table);

	Interval S;
	intRange.sup(S);

	return S;
}

Interval rho(const TaylorModelVec & tmv, const RowVector & l, const vector<Interval> & domain)
{
	int d = l.size();
	TaylorModel tmObj;

	for(int i=0; i<d; ++i)
	{
		TaylorModel tmTemp;
		Interval intTemp(l.get(i));
		tmv.tms[i].mul(tmTemp, intTemp);
		tmObj.add_assign(tmTemp);
	}

	Interval intRange;
	tmObj.intEval(intRange, domain);

	Interval S;
	intRange.sup(S);

	return S;
}

Interval rhoNormal(const TaylorModelVec & tmv, const RowVector & l, const vector<Interval> & step_exp_table)
{
	int d = l.size();
	TaylorModel tmObj;

	for(int i=0; i<d; ++i)
	{
		TaylorModel tmTemp;
		Interval intTemp(l.get(i));
		tmv.tms[i].mul(tmTemp, intTemp);
		tmObj.add_assign(tmTemp);
	}

	Interval intRange;
	tmObj.intEvalNormal(intRange, step_exp_table);

	Interval S;
	intRange.sup(S);

	return S;
}

void templatePolyhedron(Polyhedron & result, const TaylorModelVec & tmv, const vector<Interval> & domain)
{
	list<LinearConstraint>::iterator iter;
	for(iter=result.constraints.begin(); iter!=result.constraints.end(); ++iter)
	{
		iter->B = rho(tmv, iter->A, domain);
	}
}

void templatePolyhedronNormal(Polyhedron & result, const TaylorModelVec & tmv, vector<Interval> & step_exp_table)
{
	list<LinearConstraint>::iterator iter;
	for(iter=result.constraints.begin(); iter!=result.constraints.end(); ++iter)
	{
		iter->B = rhoNormal(tmv, iter->A, step_exp_table);
	}
}

int intersection_check_interval_arithmetic(const vector<PolynomialConstraint> & pcs, const vector<HornerForm> & objFuncs, const vector<Interval> & remainders, const vector<Interval> & domain, vector<bool> & bNeeded)
{
	int counter = 0;
	bNeeded.clear();

	for(int i=0; i<pcs.size(); ++i)
	{
		Interval intTemp;
		objFuncs[i].intEval(intTemp, domain);
		intTemp += remainders[i];

		if(intTemp > pcs[i].B)
		{
			// no intersection
			return -1;
		}
		else if(intTemp.smallereq(pcs[i].B))
		{
			// the flowpipe is entirely contained in the guard, domain contraction is not needed
			bNeeded.push_back(false);
			++counter;
		}
		else
		{
			bNeeded.push_back(true);
		}
	}

	return counter;
}

bool boundary_intersected_collection(const vector<PolynomialConstraint> & pcs, const vector<HornerForm> & objFuncs, const vector<Interval> & remainders, const vector<Interval> & domain, vector<bool> & boundary_intersected)
{
	boundary_intersected.clear();

	for(int i=0; i<pcs.size(); ++i)
	{
		Interval intTemp;
		objFuncs[i].intEval(intTemp, domain);
		intTemp += remainders[i];

		if(intTemp > pcs[i].B)
		{
			// no intersection
			return false;
		}
		else if(intTemp < pcs[i].B)
		{
			// do not intersect the boundary
			boundary_intersected.push_back(false);
		}
		else
		{
			boundary_intersected.push_back(true);
		}
	}

	return true;
}

int contract_interval_arithmetic(TaylorModelVec & flowpipe, vector<Interval> & domain, const vector<PolynomialConstraint> & pcs, vector<bool> & boundary_intersected)
{
	int rangeDim = flowpipe.tms.size();
	int domainDim = domain.size();

	// the Horner forms of p(T(x))
	vector<HornerForm> objHF;
	vector<Interval> remainders;

	vector<Interval> flowpipePolyRange;
	flowpipe.polyRange(flowpipePolyRange, domain);

	for(int i=0; i<pcs.size(); ++i)
	{
		TaylorModel tmTemp;
		pcs[i].hf.insert(tmTemp, flowpipe, flowpipePolyRange, domain);

		HornerForm hf;
		Interval remainder;
		tmTemp.toHornerForm(hf, remainder);
		objHF.push_back(hf);
		remainders.push_back(remainder);
	}

	// check the type of the intersections
	vector<bool> bNeeded;

	int counter = intersection_check_interval_arithmetic(pcs, objHF, remainders, domain, bNeeded);

	if(counter == -1)
	{
		return -1;	// no intersection is detected
	}
	else if(counter == pcs.size())
	{
		return 0;	// domain is not contracted
	}

	Interval intTime = domain[0];

	bool bvalid = true;
	bool bcontinue = true;

	Interval W;

	for(; bcontinue; )
	{
		vector<Interval> oldDomain = domain;

		// contract the domain
		for(int i=0; i<domainDim; ++i)
		{
			Interval newInt = domain[i];
			vector<bool> localNeeded = bNeeded;
			int localCounter = counter;

			newInt.width(W);

			// search an approximation for the lower bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<pcs.size(); ++j)
				{
					if(localNeeded[j])
					{
						vector<Interval> newDomain = domain;
						newDomain[i] = intLeft;

						Interval intTemp;
						objHF[j].intEval(intTemp, newDomain);
						intTemp += remainders[j];

						if(intTemp > pcs[j].B)
						{
							// no intersection on the left half
							newInt = intRight;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(pcs[j].B))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intLeft;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intLeft;
							newInt.width(W);

							continue;
						}
					}
				}

				if(localCounter == pcs.size())
				{
					break;
				}
			}

			// set the lower bound
			Interval Inf;
			newInt.inf(Inf);
			domain[i].setInf(Inf);

			newInt = domain[i];

			localNeeded = bNeeded;
			localCounter = counter;

			newInt.width(W);

			// search an approximation for the upper bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<pcs.size(); ++j)
				{
					if(localNeeded[j])
					{
						vector<Interval> newDomain = domain;
						newDomain[i] = intRight;

						Interval intTemp;
						objHF[j].intEval(intTemp, newDomain);
						intTemp += remainders[j];

						if(intTemp > pcs[j].B)
						{
							// no intersection on the right half
							newInt = intLeft;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(pcs[j].B))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intRight;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intRight;
							newInt.width(W);
							continue;
						}
					}
				}

				if(localCounter == pcs.size())
				{
					break;
				}
			}

			Interval Sup;
			newInt.sup(Sup);
			domain[i].setSup(Sup);	// set the upper bound

			if(!domain[i].valid())
			{
				bvalid = false;
				break;
			}
		}

		if(!bvalid)
		{
			break;
		}

		bcontinue = false;
		for(int i=0; i<domainDim; ++i)
		{
			if(oldDomain[i].widthRatio(domain[i]) <= DC_THRESHOLD_IMPROV)
			{
				bcontinue = true;
				break;
			}
		}

		if(bcontinue)
		{
			objHF.clear();
			remainders.clear();

			flowpipe.polyRange(flowpipePolyRange, domain);

			for(int i=0; i<pcs.size(); ++i)
			{
				TaylorModel tmTemp;

				pcs[i].hf.insert(tmTemp, flowpipe, flowpipePolyRange, domain);

				HornerForm hf;
				Interval remainder;
				tmTemp.toHornerForm(hf, remainder);
				objHF.push_back(hf);
				remainders.push_back(remainder);
			}
		}
	}

	if(bvalid)
	{
		boundary_intersected_collection(pcs, objHF, remainders, domain, boundary_intersected);
	}

	// normalize the contracted flowpipe
	flowpipe.normalize(domain);

	if(!bvalid)
	{
		return -1;
	}

	if(intTime != domain[0])
	{
		return 2;
	}
	else
	{
		return 1;
	}
}

void gridBox(list<vector<Interval> > & grids, const vector<Interval> & box, const int num)
{
	grids.clear();
	grids.push_back(box);

	for(int i=0; i<box.size(); ++i)
	{
		list<vector<Interval> >::iterator gridIter;
		list<vector<Interval> > newGrids;

		for(; grids.size() > 0;)
		{
			gridIter = grids.begin();

			list<Interval> queue;
			(*gridIter)[i].split(queue, num);

			list<Interval>::iterator iterComponent = queue.begin();
			for(; iterComponent != queue.end(); ++iterComponent)
			{
				vector<Interval> tmpBox = *gridIter;
				tmpBox[i] = *iterComponent;
				newGrids.push_back(tmpBox);
			}

			grids.pop_front();
		}

		grids = newGrids;
	}
}


