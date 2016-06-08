#!/bin/bash

#
# script for testing convergence rate of scheme on sample problem
# ensure diffusion is compiled with TEST_CASE = 1
#

HOSTNAME=$(hostname)
EXEC="diffusion-test"
DATE=$(date)
USER=$(whoami)

CMAKE_SRC_DIR=$(cat ../../../CMakeCache.txt | grep -E '.*_SOURCE_DIR:STATIC=.*' | grep -E -o '/.*')
GITLOG=$(cd ${CMAKE_SRC_DIR} && git log -1 | sed 's/^./# &/g' | sed 's/^$/# &/g')

repeats="5"
tau="1.0e-1"
endtime="0.5"
level="0"
FLAGS="$@\
 fem.io.macroGridFile_2d:../data/sphere.dgf\
 diffusion.repeats:${repeats}\
 diffusion.timestep:${tau}\
 diffusion.endtime:${endtime}\
 diffusion.level:${level}\
 diffusion.reducetimestepfactor:0.25\
 diffusion.D:1.0\
 diffusion.problem:diffusion\
 diffusion.streamlinediffusion:1 \
 diffusion.subtracttangentialvelocity:1 \
 fem.io.savecount:-1\
 fem.io.savestep:0.05\
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

