#!/bin/bash
# environment settings for use on Juno HPC cluster
# USAGE: . env.juno.sh <target>

# set -eu # NOTE: do not use this because it can propagate into cwl-tool subprocesses and make them break unneccessarily
arg="${1:-None}"

case $arg in
    test)
        module load singularity/3.3.0
        module load python/3.7.1
        module load cwl/cwltool
        ;;

    singularity)
        unset SINGULARITY_CACHEDIR  # TODO: why do I need to unset this when I am explicitly setting it in the Makefile??
        module load singularity/3.3.0
        ;;

    shell)
        module load singularity/3.3.0
        module load python/3.7.1
        module load cwl/cwltool
        module load jq
        ;;

    conda)
        module load singularity/3.3.0
        export PATH=${PWD}/conda/bin:${PATH}
        unset PYTHONPATH
        unset PYTHONHOME

  # *)
  #   echo "something else was called"
  #   ;;
esac
