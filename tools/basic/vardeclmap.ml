open Batteries

type var = Vardecl.var
type value = Vardecl.value
type vardecl = Vardecl.t
type t = (var, value) Map.t

let of_list (vardecls : vardecl list) : t
  =
  List.fold_left
    (fun (map : t) ((var, value) : vardecl) ->
     Map.add var value map
    )
    Map.empty
    vardecls

let print out = Map.print String.print Value.print out

let find key map =
  try
    Map.find key map
  with e ->
    let out = IO.stderr in
    begin
      String.println out "Vardeclmap Exception!";
      Printf.fprintf out "Key: %s\n" key;
      Printf.fprintf out "Map: %s\n" (IO.to_string print map);
      Printexc.print_backtrace IO.stdout;
      raise e
    end
