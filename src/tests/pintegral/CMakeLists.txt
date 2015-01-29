set(DREAL_BINARY ${DREAL_SOURCE_DIR}/../bin/dReal)
file(GLOB_RECURSE PINTEGRAL_INSTANCES ${DREAL_SOURCE_DIR}/tests/pintegral "*.smt2")
list(SORT PINTEGRAL_INSTANCES)
FOREACH(instance ${PINTEGRAL_INSTANCES})
  GET_FILENAME_COMPONENT(instance_name ${instance} NAME)
  add_test(NAME "regression_pintegral_${instance_name}"
           COMMAND "${DREAL_SOURCE_DIR}/tests/pintegral/test.sh" ${DREAL_BINARY} ${instance})
  set_tests_properties("regression_pintegral_${instance_name}" PROPERTIES LABELS "pintegral_normal")
ENDFOREACH(instance)
