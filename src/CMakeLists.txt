add_subdirectory("rg_mesh")
add_subdirectory("triangulated-surfaces")

add_executable("dune-imagedata" dune-imagedata.cc)
target_link_dune_default_libraries("dune-imagedata")
