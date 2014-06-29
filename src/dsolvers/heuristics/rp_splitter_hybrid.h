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

#include "realpaver/rp_problem.h"
#include "realpaver/rp_box.h"
#include "realpaver/rp_box_set.h"
#include "realpaver/rp_split_selector.h"
#include "realpaver/rp_split.h"
#include "dsolvers/heuristics/ode_sim_heuristic.h"

using dreal::ode_sim_heuristic;

// -------------------------------------------
// Class for domain splitting
// Enumeration of integer, bisection of reals
// Uses Value Selection Heuristic for Hybrid
// Systems
// -------------------------------------------
class rp_splitter_mixed_hybrid : public rp_splitter {
public:
  // Constructor
  rp_splitter_mixed_hybrid(rp_problem * p);

  // Destructor
  ~rp_splitter_mixed_hybrid();

  // Split the current box from bs along dimension var
  void apply(rp_box_set& bs, int var);

  void initialize(ode_sim_heuristic &ode_sim_heuristic);


private:
  // Copy protection
  rp_splitter_mixed_hybrid(const rp_splitter_mixed_hybrid& ds);
  rp_splitter_mixed_hybrid& operator=(const rp_splitter_mixed_hybrid& ds);

  ode_sim_heuristic *m_ode_sim_heuristic;
};

