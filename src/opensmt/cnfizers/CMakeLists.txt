add_library(opensmt_cnfizers Cnfizer.C Tseitin.C)
add_dependencies(opensmt_cnfizers util)