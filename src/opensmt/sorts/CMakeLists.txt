add_library(opensmt_sorts SStore.C Snode.C)
target_link_libraries(opensmt_sorts util)
add_dependencies(opensmt_sorts util)