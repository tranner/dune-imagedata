set(EXTRA_DIST  parameter *.dgf)

set(DGFFILES sphere.dgf sphere2.dgf sphere-0.3.dgf )

dune_symlink_to_source_files(FILES parameter ${DGFFILES})
