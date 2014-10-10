/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#include "Monomial.h"

double cutoff_threshold;

Monomial::Monomial()
{
}

Monomial::Monomial(const Interval & I, const vector<int> & degs):coefficient(I), degrees(degs), d(0)
{
	for(int i=0; i<degs.size(); ++i)
	{
		d += degs[i];
	}
}

Monomial::Monomial(const Monomial & monomial): coefficient(monomial.coefficient), degrees(monomial.degrees), d(monomial.d)
{
}

Monomial::Monomial(const Interval & I, const int numVars):d(0)
{
	for(int i=0; i<numVars; ++i)
	{
		degrees.push_back(0);
	}

	coefficient = I;
}

Monomial::~Monomial()
{
	degrees.clear();
}

int Monomial::degree() const
{
	return d;
}

int Monomial::dimension() const
{
	return degrees.size();
}

void Monomial::intEval(Interval & result, const vector<Interval> & domain) const
{
	result = coefficient;

	for(int i=0; i<degrees.size(); ++i)
	{
		Interval tmpI(1,1);
		for(int j=0; j<degrees[i]; ++j)
		{
			tmpI *= domain[i];
		}
		result *= tmpI;
	}
}

void Monomial::intEvalNormal(Interval & result, const vector<Interval> & step_exp_table) const
{
	Interval intZero;
	result = intZero;

	if(degrees.size() == 0)
		return;

	result = coefficient;
	result *= step_exp_table[degrees[0]];

	Interval evenInt(0,1), oddInt(-1,1);
	Interval intFactor(1);
	bool bSet = false;

	for(int i=1; i<degrees.size(); ++i)
	{
		if(degrees[i] == 0)			// degree is zero
		{
			continue;
		}
		else if(degrees[i]%2 == 0)	// degree is an even number
		{
			if(!bSet)
			{
				intFactor = evenInt;
				bSet = true;
			}
		}
		else						// degree is an odd number
		{
			intFactor = oddInt;
			break;
		}
	}

	result *= intFactor;
}

void Monomial::inv(Monomial & result) const
{
	result = *this;
	coefficient.inv(result.coefficient);
}

Monomial & Monomial::operator = (const Monomial & monomial)
{
	if(this == &monomial)
		return *this;

	coefficient = monomial.coefficient;
	degrees = monomial.degrees;
	d = monomial.d;

	return *this;
}

Monomial & Monomial::operator += (const Monomial & monomial)
{
	coefficient += monomial.coefficient;
	return *this;
}

Monomial & Monomial::operator *= (const Monomial & monomial)
{
	coefficient *= monomial.coefficient;

	for(int i=0; i<degrees.size(); ++i)
	{
		degrees[i] += monomial.degrees[i];
	}

	d += monomial.d;
	return *this;
}

const Monomial Monomial::operator + (const Monomial & monomial) const
{
	Monomial result = *this;
	result += monomial;
	return result;
}

const Monomial Monomial::operator * (const Monomial & monomial) const
{
	Monomial result = *this;
	result *= monomial;
	return result;
}

bool Monomial::isLinear(int & index) const
{
	if(d == 1)
	{
		for(int i=0; i<degrees.size(); ++i)
		{
			if(degrees[i] == 1)
			{
				index = i;
				return true;
			}
		}
	}

	return false;
}

void Monomial::dump_interval(FILE *fp, const vector<string> & varNames) const
{
	coefficient.dump(fp);

	for(int i=0; i<degrees.size()-1; i++)
	{
		if(degrees[i] != 0)
		{
			if(degrees[i] == 1)
				fprintf(fp, " * %s", varNames[i].c_str());
			else
				fprintf(fp, " * %s^%d", varNames[i].c_str(), degrees[i]);
		}
	}

	if(degrees[degrees.size()-1] != 0)
	{
		if(degrees[degrees.size()-1] == 1)
			fprintf(fp, " * %s", varNames[degrees.size()-1].c_str());
		else
			fprintf(fp, " * %s^%d", varNames[degrees.size()-1].c_str(), degrees[degrees.size()-1]);
	}
}

void Monomial::dump_constant(FILE *fp, const vector<string> & varNames) const
{
	double c = coefficient.sup();
	fprintf(fp, "(%lf)", c);

	for(int i=0; i<degrees.size()-1; i++)
	{
		if(degrees[i] != 0)
		{
			if(degrees[i] == 1)
				fprintf(fp, " * %s", varNames[i].c_str());
			else
				fprintf(fp, " * %s^%d", varNames[i].c_str(), degrees[i]);
		}
	}

	if(degrees[degrees.size()-1] != 0)
	{
		if(degrees[degrees.size()-1] == 1)
			fprintf(fp, " * %s", varNames[degrees.size()-1].c_str());
		else
			fprintf(fp, " * %s^%d", varNames[degrees.size()-1].c_str(), degrees[degrees.size()-1]);
	}
}

void Monomial::toString(string & result, const vector<string> & varNames) const
{
	string strMono;

	strMono += '(';

	string strInt;
	coefficient.toString(strInt);
	strMono += strInt;

	for(int i=0; i<degrees.size(); i++)
	{
		if(degrees[i] != 0)
		{
			if(degrees[i] == 1)
			{
				strMono += ' ';
				strMono += '*';
				strMono += ' ';
				strMono += varNames[i];
			}
			else
			{
				strMono += ' ';
				strMono += '*';
				strMono += ' ';
				strMono += varNames[i];
				strMono += '^';

				char strNum[NUM_LENGTH];
				sprintf(strNum, "%d", degrees[i]);
				string num(strNum);
				strMono += num;
			}
		}
	}

	strMono += ')';

	result = strMono;
}

bool Monomial::classInvariantOK() const
{
	int sum = 0;

	for(int i = 0; i<degrees.size(); ++i)
		sum += degrees[i];
	return (sum == d);
}

bool Monomial::cutoff(Monomial & monoRem)
{
	Interval M, intCutoff(-cutoff_threshold, cutoff_threshold);
	coefficient.midpoint(M);

	monoRem = *this;

	if(M.subseteq(intCutoff))
	{
		return false;
	}
	else
	{
		monoRem.coefficient -= M;
		coefficient = M;
		return true;
	}
}

bool Monomial::cutoff()
{
	Interval M, intCutoff(-cutoff_threshold, cutoff_threshold);
	coefficient.midpoint(M);

	if(M.subseteq(intCutoff))
	{
		return false;
	}
	else
	{
		coefficient = M;
		return true;
	}
}

bool operator == (const Monomial & a, const Monomial & b)
{
	if (a.d == b.d)
	{
		for(int i=0; i<a.degrees.size(); i++)
		{
			if(a.degrees[i] != b.degrees[i])
				return false;
		}
		return true;	// The two monomials are identical without considering the coefficients.
	}
	else
		return false;
}

bool operator < (const Monomial & a, const Monomial & b)
{
	if(a.d < b.d)
		return true;
	else if(a.d > b.d)
		return false;
	else	// a.d == b.d
	{
		for(int i=0; i<a.degrees.size(); ++i)
		{
			if(a.degrees[i] < b.degrees[i])
				return true;
			else if(a.degrees[i] > b.degrees[i])
				return false;
		}
	}
	return false;	// a == b
}

