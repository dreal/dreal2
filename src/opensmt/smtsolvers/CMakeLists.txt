add_library(opensmt_smtsolvers CoreSMTSolver.C SimpSMTSolver.C Proof.C Debug.C Theory.C TheoryCombination.C SMTConfig.C)
target_link_libraries(opensmt_smtsolvers util gflags opensmt_heuristics)
add_dependencies(opensmt_smtsolvers util)
add_dependencies(opensmt_smtsolvers GFLAGS)