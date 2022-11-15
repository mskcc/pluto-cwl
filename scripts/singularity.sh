#!/bin/bash
# test bed wrapper script around Singularity in order to develop new CWL workflow components
set -eux

SIF="${SIF:-/juno/work/ci/singularity_cache_helix_filters/pull/mskcc_helix_filters_01:21.4.1.sif}"

module load singularity/3.7.1

singularity exec \
--bind /juno \
--bind ${PWD} \
${SIF} \
bash run.sh
# bash test2.sh

# --quiet \
# --ipc \
# --pid \
# --home \
# --contain \
# --pwd /kkXAtc \
