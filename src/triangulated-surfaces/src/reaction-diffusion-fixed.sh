#!/bin/bash

HOSTNAME=`hostname`
# fem.io.macroGridFile_2d:../../rg_mesh/vtkpipe/stack_rgmesh000001_scaled.dgf \

# reactiondiffusion.level:5 \

# KochMeinhardt.D1:0.005 KochMeinhardt.D2:0.2 \
# KochMeinhardt.a0:1.0 KochMeinhardt.s0:1.25 \
# KochMeinhardt.sigma1:0.0 KochMeinhardt.sigma2:2.0e-2 \
# KochMeinhardt.mu1:1.0e-2 KochMeinhardt.mu2:0.0 \
# KochMeinhardt.rho1:1.0e-2 KochMeinhardt.rho2:2.0e-2 \
# KochMeinhardt.K:0.25 \
# KochMeinhardt.initialvariation:0.1 

# grayscott.D1:1 \
# grayscott.D2:2 \
# grayscott.F:0.0580 \
# grayscott.k:0.0630 \


# fem.io.macroGridFile_2d:../data/sphere-0.3.dgf \

# FLAGS="$@ \
# reactiondiffusion.level:5 \
# fem.io.macroGridFile_2d:../data/sphere-0.3.dgf \
# reactiondiffusion.problem:AllenCahn \
# reactiondiffusion.timestep:1.0e-4 \
# fem.io.savecount:-1 fem.io.savestep:0.025 \
# reactiondiffusion.endtime:10.0 \
# allencahn.eps:5.0e-3 \
# "

DATE=`date`
USER=`whoami`

CMAKE_SRC_DIR=$(cat ../../../CMakeCache.txt | grep -E '.*_SOURCE_DIR:STATIC=.*' | grep -E -o '/.*')
GITLOG=$(cd ${CMAKE_SRC_DIR} && git log -1 | sed 's/^./# &/g' | sed 's/^$/# &/g')

WHEN=`date +%d%m%y_%H%M%S`
FN="$0.${WHEN}"
exec  > >(tee pbs/${FN}.out) 2> >(tee pbs/${FN}.err >&2)

for f in 001 050 150 200
do
EXEC="reaction-diffusion-fixed-mesh"
FLAGS="$@ \
fem.io.macroGridFile_2d:../../rg_mesh/vtkpipe/stack_rgmesh000${f}_scaled.dgf \
reactiondiffusion.timestep:1.0e-4 \
fem.io.savecount:-1 fem.io.savestep:0.1 \
reactiondiffusion.problem:Brusselator \
reactiondiffusion.endtime:10.0
Brusselator.D1:1.0 \
Brusselator.D2:10.0 \
Brusselator.alpha:0.1 \
Brusselator.beta:0.9 \
Brusselator.gamma:200.0 \
Brusselator.r1:0.0 \
Brusselator.r2:0.0 \
fem.prefix:../output/reaction-diffusion-stationary-${f} \
"

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

