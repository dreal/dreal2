/*********************************************************************
Author: Daniel Bryce <dbryce@sift.net>
        Soonho Kong <soonhok@cs.cmu.edu>
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
#include <unordered_map>
#include "opensmt/smtsolvers/SMTConfig.h"
#include "util/scoped_vec.h"
#include "dsolvers/heuristics/ode_mode_sim.h"
#include "realpaver/realpaver.h"

namespace dreal {
class ode_sim_heuristic {
public:
    ode_sim_heuristic(){}
    void initialize(rp_propagator &,  rp_problem &);
    ~ode_sim_heuristic() {
    }

    void add_mode(SMTConfig & config, Egraph & egraph, Enode* l_int, // vector<Enode*> invs,
                  std::unordered_map<Enode*, int> & enode_to_rp_id,
                  std::unordered_map<int, Enode*> & rp_id_to_enode);
    rp_box sim(rp_box &box, int varToGet);
    ode_sim_heuristic& operator=(const ode_sim_heuristic& ds);
    void pprint_vars(ostream & out, rp_box b, bool exact) const;
    void display_interval(ostream & out, rp_interval i, int digits, bool exact) const;
private:
    SMTConfig m_config;
    rp_propagator *m_propag;
    vector<Enode*> m_integral_lits;
    map< int, ode_mode_sim* > m_mode_sims;
    rp_problem *m_rp_problem;
    std::unordered_map<int, Enode*> * m_rp_id_to_enode;
    std::unordered_map<Enode *, int> * m_enode_to_rp_id;
};
}
