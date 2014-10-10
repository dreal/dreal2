/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#ifndef CONTINUOUS_H_
#define CONTINUOUS_H_

#include "TaylorModel.h"
#include "Geometry.h"

class Flowpipe					// A flowpipe is represented by a composition of two Taylor models. The left Taylor model is the preconditioning part.
{
private:
	TaylorModelVec tmvPre;		// preconditioning Taylor model
	TaylorModelVec tmv;
private:
	vector<Interval> domain;	// domain of TMV_right, the first variable is t
public:
	Flowpipe();
	Flowpipe(const TaylorModelVec & tmvPre_input, const TaylorModelVec & tmv_input, const vector<Interval> & domain_input);
	Flowpipe(const vector<Interval> & box, const Interval & I);								// represent a box
	Flowpipe(const TaylorModelVec & tmv_input, const vector<Interval> & domain_input);		// construct a flowpipe from a Taylor model
	Flowpipe(const Flowpipe & flowpipe);
	~Flowpipe();

	void clear();
	void dump(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames) const;
	void dump_normal(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames, vector<Interval> & step_exp_table) const;
	void composition(TaylorModelVec & result) const;	// apply the preconditioning part to the Taylor model
	void composition_normal(TaylorModelVec & result, const vector<Interval> & step_exp_table) const;

	void intEval(vector<Interval> & result) const;
	void intEvalNormal(vector<Interval> & result, const vector<Interval> & step_exp_table) const;

	void normalize();

	// fast integration scheme for low-degree ODEs
	// fixed step sizes and orders
	bool advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;
	bool advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;

	// adaptive step sizes and fixed orders
	bool advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;
	bool advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;

	// adaptive orders and fixed step sizes
	bool advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;
	bool advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;

	// integration scheme for high-degree ODEs
	// fixed step sizes and orders
	bool advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;
	bool advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;

	// adaptive step sizes and fixed orders
	bool advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;
	bool advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;

	// adaptive orders and fixed step sizes
	bool advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;
	bool advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const int localMaxOrder, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<Interval> & uncertainties) const;



	// integration scheme for non-polynomial ODEs (using Taylor approximations)
	// fixed step sizes and orders
	bool advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const;
	bool advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const;

	// adaptive step sizes and fixed orders
	bool advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const;
	bool advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const;

	// adaptive orders and fixed step sizes
	bool advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const;
	bool advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const int localMaxOrder, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<Interval> & uncertainties, const vector<Interval> & uncertainty_centers) const;

	Flowpipe & operator = (const Flowpipe & flowpipe);

	friend class ContinuousSystem;
	friend class ContinuousReachability;
	friend class HybridSystem;
	friend class HybridReachability;
};

class ContinuousSystem
{
private:
	TaylorModelVec tmvOde;
	vector<HornerForm> hfOde;		// a Horner form of the ode
	Flowpipe initialSet;			// the initial set
	vector<Interval> uncertainties;
	vector<Interval> uncertainty_centers;
	vector<string> strOde;
public:
	ContinuousSystem();
	ContinuousSystem(const TaylorModelVec & ode_input, const vector<Interval> & uncertainties_input, const Flowpipe & initialSet_input);
	ContinuousSystem(const vector<string> & strOde_input, const vector<Interval> & uncertainties_input, const Flowpipe & initialSet_input);
	ContinuousSystem(const ContinuousSystem & system);
	~ContinuousSystem();

	// for low-degree ODEs
	// fixed step sizes and orders
	void reach_low_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;
	void reach_low_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;

	// adaptive step sizes and fixed orders
	void reach_low_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;
	void reach_low_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;

	// adaptive orders and fixed step sizes
	void reach_low_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;
	void reach_low_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;

	// for high-degree ODEs
	// fixed step sizes and orders
	void reach_high_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;
	void reach_high_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;

	// adaptive step sizes and fixed orders
	void reach_high_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;
	void reach_high_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;

	// adaptive orders and fixed step sizes
	void reach_high_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;
	void reach_high_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;

	// for non-polynomial ODEs (using Taylor approximations)
	// fixed step sizes and orders
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;

	// adaptive step sizes and fixed orders
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;

	// adaptive orders and fixed step sizes
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames) const;

	ContinuousSystem & operator = (const ContinuousSystem & system);
};

class ContinuousReachability		// The reachability analysis of continuous systems
{
public:
	ContinuousSystem system;		// the continuous system
	double step;					// the step size used in the reachability analysis
	double time;					// the time horizon for the reachability analysis
	int precondition;				// the preconditioning technique
	vector<int> outputAxes;			// the output axes
	int plotSetting;
	int plotFormat;
	int numSections;				// the number of sections in each dimension

	int orderType;
	bool bAdaptiveSteps;
	bool bAdaptiveOrders;

	vector<Interval> estimation;	// the remainder estimation for varying time step
	double miniStep;				// the minimum step size
	vector<int> orders;				// the order(s)
	vector<int> maxOrders;			// the maximum orders
	int globalMaxOrder;

	bool bPrint;
	bool bSafetyChecking;

	int integrationScheme;

	list<Flowpipe> flowpipes;
	list<TaylorModelVec> flowpipesCompo;
	list<vector<Interval> > domains;

	vector<PolynomialConstraint> unsafeSet;

	map<string,int> stateVarTab;
	vector<string> stateVarNames;

	map<string,int> tmVarTab;
	vector<string> tmVarNames;

	char outputFileName[NAME_SIZE];
public:
	ContinuousReachability();
	~ContinuousReachability();

	void dump(FILE *fp) const;

	void run();
	void composition();
	int safetyChecking() const;
	unsigned long numOfFlowpipes() const;

	void dump_potential_counterexample(FILE *fp, const list<TaylorModelVec> & flowpipes, const list<vector<Interval> > & domains, const list<Interval> & globalTimes) const;

	void plot_2D() const;

	void plot_2D_GNUPLOT(FILE *fp) const;
	void plot_2D_interval_GNUPLOT(FILE *fp) const;
	void plot_2D_octagon_GNUPLOT(FILE *fp) const;
	void plot_2D_grid_GNUPLOT(FILE *fp) const;

	void plot_2D_MATLAB(FILE *fp) const;
	void plot_2D_interval_MATLAB(FILE *fp) const;
	void plot_2D_octagon_MATLAB(FILE *fp) const;
	void plot_2D_grid_MATLAB(FILE *fp) const;

	bool declareStateVar(const string & vName);
	int getIDForStateVar(const string & vName) const;
	bool getStateVarName(string & vName, const int id) const;

	bool declareTMVar(const string & vName);
	int getIDForTMVar(const string & vName) const;
	bool getTMVarName(string & vName, const int id) const;
};

void computeTaylorExpansion(TaylorModelVec & result, const TaylorModelVec & first_order_deriv, const TaylorModelVec & ode, const int order);
void computeTaylorExpansion(TaylorModelVec & result, const TaylorModelVec & first_order_deriv, const TaylorModelVec & ode, const vector<int> & orders);

void construct_step_exp_table(vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const int order);
void construct_step_exp_table(vector<Interval> & step_exp_table, const Interval & step, const int order);

void preconditionQR(Matrix & result, const TaylorModelVec & tmv, const int rangeDim, const int domainDim);

Interval rho(const TaylorModelVec & tmv, const vector<Interval> & l, const vector<Interval> & domain);
Interval rhoNormal(const TaylorModelVec & tmv, const vector<Interval> & l, const vector<Interval> & step_end_exp_table);

Interval rho(const TaylorModelVec & tmv, const RowVector & l, const vector<Interval> & domain);
Interval rhoNormal(const TaylorModelVec & tmv, const RowVector & l, const vector<Interval> & step_end_exp_table);

void templatePolyhedron(Polyhedron & result, const TaylorModelVec & tmv, const vector<Interval> & domain);
void templatePolyhedronNormal(Polyhedron & result, const TaylorModelVec & tmv, vector<Interval> & step_end_exp_table);

int intersection_check_interval_arithmetic(const list<PolynomialConstraint> & pcs, const list<HornerForm> & objFuncs, const list<Interval> & remainders, const vector<Interval> & domain, list<bool> & bNeeded);
bool boundary_intersected_collection(const vector<PolynomialConstraint> & pcs, const vector<HornerForm> & objFuncs, const vector<Interval> & remainders, const vector<Interval> & domain, vector<bool> & boundary_intersected);

// domain contraction by using interval arithmetic
int contract_interval_arithmetic(TaylorModelVec & flowpipe, vector<Interval> & domain, const vector<PolynomialConstraint> & pcs, vector<bool> & boundary_intersected);

void gridBox(list<vector<Interval> > & grids, const vector<Interval> & box, const int num);

#endif /* CONTINUOUS_H_ */
