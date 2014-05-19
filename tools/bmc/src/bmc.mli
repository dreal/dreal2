(*
 * Soonho Kong soonhok@cs.cmu.edu
 * Wei Chen weichen1@andrew.cmu.edu
 *)
open Type
open Type.Hybrid
open Type.Basic
open Type.Mode
open Type.Jump
(*open Heuristic
open Heuristic.Costmap*)
open Batteries
open IO
open Smt2_cmd
open Smt2

exception SMTException of string

(** a list of annoted flow ode **)
type flows_annot = (int * ode list)  (** step, mode, ode **)

(** compile a Hybrid automata into SMT formula **)
val compile : Hybrid.t -> int -> int list option -> Smt2.t
val pathgen : Hybrid.t -> int -> (int list) list
(* val heuristicgen : Hybrid.t -> int -> Costmap.t *)
