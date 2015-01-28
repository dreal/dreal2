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

#include "dsolvers/heuristics/rp_splitter_time.h"
#include <algorithm>
#include <set>
#include <map>
#include <iostream>


// -------------------------------------------
// Class for domain splitting
// Enumeration of integer, bisection of reals
// Uses Hybrid System value selection heuristic
// -------------------------------------------
rp_splitter_time::rp_splitter_time(rp_problem * p):
  rp_splitter(p)
{}


void rp_splitter_time::initialize(std::set<int> *tvars, double nra_precision){
  time_vars = tvars;
  precision = nra_precision/100.0;
}

rp_splitter_time::~rp_splitter_time(){
  delete time_vars;
  for (std::map<rp_box*, std::set<int>*>::iterator i = quick_split_vars.begin(); i != quick_split_vars.end(); i++){
    delete (*i).second;
  }
}

bool rp_splitter_time::is_time_variable(int var){
  // std::cout << "TIME var?" << var << std::endl;
  // for (std::set<int>::iterator i = time_vars->begin(); i != time_vars->end(); i++)
  //   std::cout << *i << std::endl;
  return (time_vars->find(var) != time_vars->end());
}

bool rp_splitter_time::did_quick_split(rp_box* box, int var){
  return (quick_split_vars[box]->find(var) != quick_split_vars[box]->end());
}

void rp_splitter_time::apply(rp_box_set& bs, int var) {
  rp_interval i1, i2;
  rp_box b1 = bs.remove_insert();
  this->observe(b1, var);
  rp_box b2 = bs.insert(b1);

  if (quick_split_vars[&b1] == NULL){
    quick_split_vars[&b1] = new std::set<int>();
  }
  quick_split_vars[&b2] = new std::set<int>(*quick_split_vars[&b1]);

  if (is_time_variable(var) && !did_quick_split(&b1, var)){
    if (this->real_hole(rp_box_elem(b1, var),
                        rp_variable_domain(rp_problem_var(*_problem, var)),
                        i1, i2)){
      rp_interval_copy(rp_box_elem(b1, var), i1);
      rp_interval_copy(rp_box_elem(b2, var), i2);
    } else {
      // std::cout << "TIME SPLIT" << std::endl;
      rp_interval &vi = rp_box_elem(b1, var);
      double split =  std::min(rp_binf(vi) + precision, rp_split_point(rp_binf(vi), rp_bsup(vi), 1000, 1));

      // Real variable: [a,b] --> [center,b] and [a,center]
      rp_binf(rp_box_elem(b1, var)) =
        rp_bsup(rp_box_elem(b2, var)) =
        split;

      quick_split_vars[&b1]->insert(var);
      quick_split_vars[&b2]->insert(var);
        // rp_interval_midpoint(rp_box_elem(b1,var));
    }
  } else if (rp_variable_integer(rp_problem_var(*_problem, var))){
    if (this->integer_hole(rp_box_elem(b1, var),
                           rp_variable_domain(rp_problem_var(*_problem, var)),
                           i1, i2)) {
      rp_interval_copy(rp_box_elem(b1, var), i1);
      rp_interval_copy(rp_box_elem(b2, var), i2);
    } else { // no hole found
      // Integer variable: [a,b] --> [a+1,b] and [a,a]
      ++rp_binf(rp_box_elem(b1, var));
      rp_bsup(rp_box_elem(b2, var)) = rp_binf(rp_box_elem(b2, var));
    }
  } else {
    if (this->real_hole(rp_box_elem(b1, var),
                        rp_variable_domain(rp_problem_var(*_problem, var)),
                        i1, i2)) {
      rp_interval_copy(rp_box_elem(b1, var), i1);
      rp_interval_copy(rp_box_elem(b2, var), i2);
    } else {
      // Real variable: [a,b] --> [center,b] and [a,center]
      rp_binf(rp_box_elem(b1, var)) =
        rp_bsup(rp_box_elem(b2, var)) =
          rp_interval_midpoint(rp_box_elem(b1, var));
    }
  }
}

rp_splitter_time::rp_splitter_time(const rp_splitter_time& ds): rp_splitter(ds)
{}

rp_splitter_time&
rp_splitter_time::operator=(const rp_splitter_time& /*ds*/){
  return( *this );
}
