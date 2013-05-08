type logic = QF_NRA
             | QF_NRA_ODE

type formula = Basic.formula

type t = SetLogic of logic
         | SetInfo of string * string
         | DeclareFun of string
         | DeclareConst of string
         | Assert of formula
         | CheckSAT
         | Exit

let logic_qfnra =
  SetLogic QF_NRA

let make_lb (name : string) (v : float)
    = Assert (Basic.Le (Basic.Num v,  Basic.Var name))
let make_ub (name : string) (v : float)
    = Assert (Basic.Le (Basic.Var name, Basic.Num v ))

let set_precision (p : float) : t =
  SetInfo (":precision", string_of_float p)

let print_logic out =
  function
  | QF_NRA -> BatString.print out "QF_NRA"
  | QF_NRA_ODE -> BatString.print out "QF_NRA_ODE"

let print out =
  function
  | SetLogic l ->
    begin
      BatString.print   out "(set-logic ";
      print_logic out l;
      BatString.print out ")";
    end
  | SetInfo (key, value) ->
    begin
      BatString.print out "(set-info ";
      BatString.print out key;
      BatString.print out " ";
      BatString.print out value;
      BatString.print out ")";
    end
  | DeclareConst v ->
    begin
      BatString.print   out "(declare-const ";
      BatString.print   out v;
      BatString.print out " Real)";
    end
  | DeclareFun v ->
    begin
      BatString.print   out "(declare-fun ";
      BatString.print   out v;
      BatString.print out " () Real)";
    end
  | Assert f ->
    begin
      BatString.print out "(assert ";
      Basic.print_formula out f;
      BatString.print out ")";
    end
  | CheckSAT ->
    BatString.print out "(check-sat)"
  | Exit ->
    BatString.print out "(exit)"