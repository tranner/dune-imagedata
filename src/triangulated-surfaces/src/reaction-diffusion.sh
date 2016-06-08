#!/bin/bash

HOSTNAME=`hostname`
DATE=`date`
USER=`whoami`

CMAKE_SRC_DIR=$(cat ../../../CMakeCache.txt | grep -E '.*_SOURCE_DIR:STATIC=.*' | grep -E -o '/.*')
GITLOG=$(cd ${CMAKE_SRC_DIR} && git log -1 | sed 's/^./# &/g' | sed 's/^$/# &/g')

WHEN=`date +%d%m%y_%H%M%S`
FN="$0.${WHEN}"
exec  > >(tee pbs/${FN}.out) 2> >(tee pbs/${FN}.err >&2)

# for initialfile in 0
# do

# EXEC="reaction-diffusion"
# FLAGS="$@ \
# reactiondiffusion.timestep:1.0e-4 \
# fem.io.savecount:-1 fem.io.savestep:0.1 \
# reactiondiffusion.problem:Brusselator \
# reactiondiffusion.endtime:100.0 \
# Brusselator.D1:1.0 \
# Brusselator.D2:10.0 \
# Brusselator.alpha:0.1 \
# Brusselator.beta:0.9 \
# Brusselator.gamma:200.0 \
# Brusselator.r1:0.0 \
# Brusselator.r2:0.0 \
# fem.prefix:../output/reaction-diffusion-evolving-${initialfile}-normal-velocity \
# imagedata.filetimestep:1.0 \
# imagedata.initialfile:${initialfile} \
# "

# echo "# ------------------------------------------------------------"
# echo "# running on ${HOSTNAME}"
# echo "# at ${DATE}"
# echo "# by ${USER}"
# echo "# calling ./${EXEC} ${FLAGS}"
# echo "# cmake_src_dir ${CMAKE_SRC_DIR}"
# echo "# most recent git log "
# echo "${GITLOG}"
# echo "# ------------------------------------------------------------"

# ./${EXEC} ${FLAGS}

# ENDDATE=`date`
# echo "# ------------------------------------------------------------"
# echo "# finised at ${ENDDATE}"
# echo "# ------------------------------------------------------------"

# done

for initialfile in 0
do

EXEC="reaction-diffusion"
FLAGS="$@ \
reactiondiffusion.timestep:1.0e-4 \
fem.io.savecount:-1 fem.io.savestep:0.1 \
reactiondiffusion.problem:Brusselator \
reactiondiffusion.endtime:100.0 \
Brusselator.D1:1.0 \
Brusselator.D2:10.0 \
Brusselator.alpha:0.1 \
Brusselator.beta:0.9 \
Brusselator.gamma:200.0 \
Brusselator.r1:0.0 \
Brusselator.r2:0.0 \
fem.prefix:../output/reaction-diffusion-evolving-${initialfile}-vertex-velocity-2 \
imagedata.filetimestep:1.0 \
imagedata.initialfile:${initialfile} \
diffusion.subtracttangentialvelocity:0 \
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

