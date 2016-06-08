#!/bin/bash

#
# script for testing convergence rate of scheme on sample problem
# ensure diffusion is compiled with TEST_CASE = 1
#

HOSTNAME=$(hostname)
EXEC="reaction-diffusion-test-fixed-mesh"
DATE=$(date)
USER=$(whoami)

CMAKE_SRC_DIR=$(cat ../../../CMakeCache.txt | grep -E '.*_SOURCE_DIR:STATIC=.*' | grep -E -o '/.*')
GITLOG=$(cd ${CMAKE_SRC_DIR} && git log -1 | sed 's/^./# &/g' | sed 's/^$/# &/g')

repeats="5"
tau="0.5"
endtime="0.5"
level="0"
FLAGS="$@\
 reactiondiffusion.repeats:${repeats} \
 reactiondiffusion.timestep:${tau} \
 reactiondiffusion.endtime:${endtime} \
 reactiondiffusion.level:${level} \
 reactiondiffusion.reducetimestepfactor:0.25 \
 fem.io.savecount:-1 \
 fem.io.savestep:-1.0 \
fem.io.macroGridFile_2d:../data/sphere.dgf \
KochMeinhardt.D1:0.005 KochMeinhardt.D2:0.2 \
KochMeinhardt.sigma1:0.0 KochMeinhardt.sigma2:2.0e-2 \
KochMeinhardt.mu1:0.0 KochMeinhardt.mu2:0.0 \
KochMeinhardt.mu1:1.0e-2 KochMeinhardt.mu2:0.0 \
KochMeinhardt.rho1:1.0e-2 KochMeinhardt.rho2:2.0e-2 \
KochMeinhardt.K:0.25 \
"





WHEN=`date +%d%m%y_%H%M%S`
FN="$0.${WHEN}"
exec  > >(tee pbs/${FN}.out) 2> >(tee pbs/${FN}.err >&2)

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

