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

#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include "dsolvers/heuristics/ode_sim_heuristic.h"
#include "opensmt/egraph/Egraph.h"
#include "util/logging.h"

using std::string;
using std::ifstream;
using std::unordered_set;
using std::ios;
using std::sort;

namespace dreal{
  void ode_sim_heuristic::initialize(rp_propagator &propag) {
    m_propag = &propag;
  }

  void ode_sim_heuristic::add_mode(SMTConfig & config, Egraph & egraph, Enode* l_int, // vector<Enode*> invs,
                                   std::unordered_map<Enode*, int> & enode_to_rp_id){
    DREAL_LOG_DEBUG << "ode_sim_heuristic::add_mode " << l_int;

    Enode* m_time = l_int->getCdr()->getCdr()->getCdr()->getCar();
    string time_str = m_time->getCar()->getName();                       // i.e. "time_1"
    int m_step = stoi(time_str.substr(time_str.find_last_of("_") + 1));      // i.e. 1
    m_mode_sims[m_step] = new ode_mode_sim(config, egraph, l_int, // invs,
                                           enode_to_rp_id);
  }

  rp_box ode_sim_heuristic::sim(rp_box &box){
    DREAL_LOG_DEBUG << "ode_sim_heuristic::sim() with modes " << m_mode_sims.size();


    rp_box & sim_box(box);


    for (int i = 0; i < static_cast<int>(m_mode_sims.size()); i++){
      DREAL_LOG_DEBUG << "ode_sim_heuristic::sim step " << i;

      vector<pair<interval, DVector>> bucket;
      ode_mode_sim & current_mode = *m_mode_sims[i];

      current_mode.update(sim_box);
      current_mode.compute_forward(bucket);

      for (auto p : bucket){
        rp_box temp_box(sim_box);
        current_mode.update_box(temp_box, p.second, p.first.rightBound());
        if (m_propag->apply(temp_box)){
          // guard is satisfied
          sim_box = temp_box;
          break;
        }
      }
    }
    return sim_box;
  }

ode_sim_heuristic&
ode_sim_heuristic::operator=(const ode_sim_heuristic& /*ds*/) {
  return( *this );
}
}
