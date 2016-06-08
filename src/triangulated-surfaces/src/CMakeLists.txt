set(GRIDTYPE  ALBERTAGRID)
set(GRIDDIM  2)
set(WORLDDIM  3)
set(POLORDER  1)
set(WANT_ISTL  1)
set(TEST_CASE 0)
set(FIXED_MESH 0)

set(CLEANFILES  *.log)

set(check_PROGRAMS  diffusion reaction-diffusion )

add_executable( diffusion diffusion.cc)
target_link_dune_default_libraries( diffusion )
add_dune_alberta_flags( diffusion GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM} )
target_compile_definitions( diffusion PRIVATE ${GRIDTYPE};GRIDDIM=${GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};TEST_CASE=${TEST_CASE};FIXED_MESH=${FIXED_MESH} )

add_executable( diffusion-fixed-mesh diffusion.cc)
target_link_dune_default_libraries( diffusion-fixed-mesh )
add_dune_alberta_flags( diffusion-fixed-mesh GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM} )
target_compile_definitions( diffusion-fixed-mesh PRIVATE ${GRIDTYPE};GRIDDIM=${GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};TEST_CASE=${TEST_CASE};FIXED_MESH=1 )

add_executable( diffusion-test diffusion.cc)
target_link_dune_default_libraries( diffusion-test )
add_dune_alberta_flags( diffusion-test GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM} )
target_compile_definitions( diffusion-test PRIVATE ${GRIDTYPE};GRIDDIM=${GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};TEST_CASE=1;FIXED_MESH=0 )

add_executable( diffusion-test-fixed-mesh diffusion.cc)
target_link_dune_default_libraries( diffusion-test-fixed-mesh )
add_dune_alberta_flags( diffusion-test-fixed-mesh GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM} )
target_compile_definitions( diffusion-test-fixed-mesh PRIVATE ${GRIDTYPE};GRIDDIM=${GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};TEST_CASE=1;FIXED_MESH=1 )

add_executable( reaction-diffusion reaction-diffusion.cc)
target_link_dune_default_libraries( reaction-diffusion )
add_dune_alberta_flags( reaction-diffusion GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM} )
target_compile_definitions( reaction-diffusion PRIVATE ${GRIDTYPE};GRIDDIM=${GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER} )

add_executable( reaction-diffusion-fixed-mesh reaction-diffusion.cc)
target_link_dune_default_libraries( reaction-diffusion-fixed-mesh )
add_dune_alberta_flags( reaction-diffusion-fixed-mesh GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM} )
target_compile_definitions( reaction-diffusion-fixed-mesh PRIVATE ${GRIDTYPE};GRIDDIM=${GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};FIXED_MESH=1 )

add_executable( reaction-diffusion-test-fixed-mesh reaction-diffusion.cc)
target_link_dune_default_libraries( reaction-diffusion-test-fixed-mesh )
add_dune_alberta_flags( reaction-diffusion-test-fixed-mesh GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM} )
target_compile_definitions( reaction-diffusion-test-fixed-mesh PRIVATE ${GRIDTYPE};GRIDDIM=${GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};FIXED_MESH=1;TEST_CASE=1 )

dune_symlink_to_source_files(FILES diffusion-frap.sh diffusion-ligand.sh reaction-diffusion.sh test.sh reaction-diffusion-test.sh reaction-diffusion-fixed.sh)
