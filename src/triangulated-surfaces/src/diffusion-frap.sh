#!/bin/bash

WHEN=`date +%d%m%y_%H%M%S`
FN="$0.${WHEN}"
exec  > >(tee pbs/${FN}.out) 2> >(tee pbs/${FN}.err >&2)

HOSTNAME=`hostname`
DATE=`date`
USER=`whoami`

CMAKE_SRC_DIR=$(cat ../../../CMakeCache.txt | grep -E '.*_SOURCE_DIR:STATIC=.*' | grep -E -o '/.*')
GITLOG=$(cd ${CMAKE_SRC_DIR} && git log -1 | sed 's/^./# &/g' | sed 's/^$/# &/g')

EXEC="diffusion"

for initialfile in 0 25 77
do
FLAGS="$@ fem.io.savecount:1 \
fem.prefix:../output/diffusion-frap-${initialfile}-normal-velocity \
 diffusion.endtime:0.0020 \
imagedata.initialfile:${initialfile}"

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

for initialfile in 0 25 77
do
FLAGS="$@ fem.io.savecount:1 \
fem.prefix:../output/diffusion-frap-${initialfile}-vertex-velocity \
diffusion.endtime:0.0020 \
imagedata.initialfile:${initialfile} \
diffusion.subtracttangentialvelocity:0"

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


