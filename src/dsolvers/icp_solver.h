/*********************************************************************
Author: Soonho Kong <soonhok@cs.cmu.edu>
        Sicun Gao <sicung@cs.cmu.edu>
        Edmund Clarke <emc@cs.cmu.edu>

dReal -- Copyright (C) 2013 - 2014, Soonho Kong, Sicun Gao, and Edmund Clarke

dReal is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

dReal is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with dReal. If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#pragma once
#include <fstream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "config.h"
#include "util/scoped_env.h"
#include "util/scoped_vec.h"
#include "opensmt/egraph/Egraph.h"
#include "opensmt/egraph/Enode.h"
#include "opensmt/smtsolvers/SMTConfig.h"
#include "realpaver/realpaver.h"
#include "dsolvers/heuristics/ode_sim_heuristic.h"

#ifdef ODE_ENABLED
#include "dsolvers/ode_solver.h"
#endif

namespace dreal {
class icp_solver {
public:
    icp_solver(SMTConfig & c, Egraph & e, SStore & t, scoped_vec const & stack, scoped_env & env, bool complete_check);
    ~icp_solver();
    bool prop(); // only propagate
    bool solve();
#ifdef ODE_ENABLED
    void        print_json(ostream& out);
#endif
    double      constraint_width(const rp_constraint * c, rp_box const b) const;
    bool        is_box_within_delta(rp_box b);
    int         get_var_split_delta(rp_box b);
    int         get_var_split_delta_hybrid(rp_box b);
    int         get_max_time_index(rp_constraint c);
    bool        delta_test() const { return m_config.nra_delta_test; }
    int         nsplit() const { return m_nsplit; }
    std::vector<Enode *> get_explanation();

private:
    icp_solver(const icp_solver& s);
    icp_solver& operator=(const icp_solver& s);
    rp_problem  create_rp_problem();
    bool        prop_with_ODE(); // propagate with ODE (only in complete check)
    rp_box      compute_next(); // computation of the next solution
    void        output_problem() const;
    void        display_box(ostream& out, rp_box b, int digits, int mode, bool exact) const;
    void        display_interval(ostream & out, rp_interval i, int digits, bool exact) const;
    void        pprint_vars(ostream & out, rp_problem p, rp_box b, bool exact) const;
    void        pprint_lits(ostream & out, rp_problem p, rp_box b) const;

#ifdef ODE_ENABLED
    void        create_ode_solvers();
    void        create_ode_sim_solvers();
    void        callODESolver(ode_solver * odeSolver, bool forward, ode_solver::ODE_result & result);
    void        print_ODE_trajectory(ostream& out) const;
#endif

    SMTConfig &                      m_config;
    Egraph &                         m_egraph;
    SStore &                         m_sstore;
    rp_problem                       m_problem; /* problem to be solved */
    rp_propagator                    m_propag; /* reduction algorithm using propagation */
    rp_box_stack                     m_boxes; /* the set of boxes during search */
    rp_selector *                    m_vselect; /* selection of variable to be split */
    rp_splitter *                    m_dsplit; /* split function of variable domain */
    int                              m_nsplit; /* number of split steps */
    double                           m_improve; /* improvement factor of iterative methods */
    vector<rp_constraint *>          m_rp_constraints;
    scoped_vec const &               m_stack;
    scoped_env &                     m_env;
    bool                             m_ODEresult;
    std::unordered_map<Enode *, int> m_enode_to_rp_id;
    std::unordered_map<Enode *, rp_constraint> m_enode_to_rp_ctr;
    std::unordered_map<int, Enode *> m_rp_id_to_enode;
#ifdef ODE_ENABLED
    std::vector<std::unordered_map<std::string, Enode *>> connect_table;
    std::vector<ode_solver *>        m_ode_solvers;
    ode_sim_heuristic                m_ode_sim_heuristic;
#endif
    bool                             m_complete_check;
    int                              m_num_delta_checks;
};
}
