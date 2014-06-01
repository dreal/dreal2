open Batteries

type formula = Basic.formula
type id = Id.t
type t = {guard  : formula;
          target: Id.t;
          change : formula}

let make (g, t, c) = {guard = g; target = t; change = c}

let guard  {guard = g; target = t; change = c} = g
let target {guard = g; target = t; change = c} = t
let change {guard = g; target = t; change = c} = c

let print out {guard  = g;
               target = t;
               change = c}
  =
  Printf.fprintf out
                 "(%s, %s, %s)"
                 (IO.to_string Basic.print_formula g)
                 (IO.to_string Id.print t)
                 (IO.to_string Basic.print_formula c)
