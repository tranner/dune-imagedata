#!/bin/bash

HOSTNAME=`hostname`
DATE=`date`
USER=`whoami`

CMAKE_SRC_DIR=$(cat ../../../CMakeCache.txt | grep -E '.*_SOURCE_DIR:STATIC=.*' | grep -E -o '/.*')
GITLOG=$(cd ${CMAKE_SRC_DIR} && git log -1 | sed 's/^./# &/g' | sed 's/^$/# &/g')

WHEN=`date +%d%m%y_%H%M%S`
FN="$0.${WHEN}"
exec  > >(tee pbs/${FN}.out) 2> >(tee pbs/${FN}.err >&2)

EXEC="diffusion-fixed-mesh"

for extra in "diffusion.ligand.cx:0 diffusion.ligand.cy:0 diffusion.ligand.cz:0.3 fem.prefix:../output/diffusion-ligand/p0 fem.io.macroGridFile_2d:../../rg_mesh/vtkpipe/stack_rgmesh000001_scaled.dgf" \
    "diffusion.ligand.cx:0 diffusion.ligand.cy:0 diffusion.ligand.cz:-0.3  fem.prefix:../output/diffusion-ligand/p1 fem.io.macroGridFile_2d:../../rg_mesh/vtkpipe/stack_rgmesh000001_scaled.dgf" \
    "diffusion.ligand.cx:0.25 diffusion.ligand.cy:0.25 diffusion.ligand.cz:0  fem.prefix:../output/diffusion-ligand/p2 fem.io.macroGridFile_2d:../../rg_mesh/vtkpipe/stack_rgmesh000001_scaled.dgf" \
    "diffusion.ligand.cx:0 diffusion.ligand.cy:0 diffusion.ligand.cz:0.1 fem.prefix:../output/diffusion-ligand/p3 fem.io.macroGridFile_2d:../../rg_mesh/vtkpipe/stack_rgmesh000150_scaled.dgf" \
    "diffusion.ligand.cx:0 diffusion.ligand.cy:0 diffusion.ligand.cz:-0.3  fem.prefix:../output/diffusion-ligand/p4 fem.io.macroGridFile_2d:../../rg_mesh/vtkpipe/stack_rgmesh000150_scaled.dgf" \
    "diffusion.ligand.cx:0.25 diffusion.ligand.cy:0.25 diffusion.ligand.cz:0  fem.prefix:../output/diffusion-ligand/p5 fem.io.macroGridFile_2d:../../rg_mesh/vtkpipe/stack_rgmesh000150_scaled.dgf"
do
    FLAGS="$@ fem.io.savestep:2.5e-4 diffusion.endtime:0.01 diffusion.problem:ligand \
diffusion.ligand.k:0 ${extra}"

    echo "# ------------------------------------------------------------"
    echo "# running on ${HOSTNAME}"
    echo "# at ${DATE}"
    echo "# by ${USER}"
    echo "# calling ./${EXEC} ${FLAGS}"
    echo "# cmake_src_dir ${CMAKE_SRC_DIR}"
    echo "# most recent git log "
    echo "${GITLOG}"
    echo "# ------------------------------------------------------------"

    ./${EXEC} ${FLAGS}

    ENDDATE=`date`
    echo "# ------------------------------------------------------------"
    echo "# finised at ${ENDDATE}"
    echo "# ------------------------------------------------------------"

done

