#!/bin/bash
# environment settings for use on Juno HPC cluster
# USAGE: . ./env.juno.sh <target>

# set -eu # NOTE: do not use this because it can propagate into cwl-tool subprocesses and make them break unneccessarily
arg="${1:-None}"

# load global 'module' command if its present
[ -e /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh || :

case $arg in
    # TODO: clean this up! Some of these are deprecated

    # for running test cases in dev with cwltool
    test)
        module load singularity/3.7.1
        module load python/3.7.1
        module load cwl/cwltool
        ;;

    # for running large end-to-end tests with Jenkins
    integration_test)
        export INTEGRATION_TESTS=True
        export USE_LSF=True
        export CWL_ENGINE=toil
        export PRINT_COMMAND="True"
        . /juno/work/ci/jenkins/pluto-cwl/toil-settings.sh
        module load singularity/3.7.1
        export PATH=/juno/work/ci/jenkins/pluto-cwl/pluto-cwl/conda/bin:${PATH}
        export PATH="$(dirname $(which singularity))":${PATH}
        unset PYTHONPATH
        unset PYTHONHOME
        ;;

    # for running all test cases including the large integration tests
    # TODO: remove this one!
    test-full)
        export LARGE_TESTS=True
        module load singularity/3.7.1
        module load python/3.7.1
        module load cwl/cwltool
        ;;

    # for running CWL's with Singularity
    # TODO: remove this one!
    singularity)
        unset SINGULARITY_CACHEDIR  # TODO: why do I need to unset this when I am explicitly setting it in the Makefile??
        module load singularity/3.7.1
        ;;

    # for interactive shell session
    # TODO: Do we still need this one??
    shell)
        module load singularity/3.7.1
        module load python/3.7.1
        module load cwl/cwltool
        module load jq
        ;;

    # for using conda install dependencies
    conda)
        export PATH=${PWD}/conda/bin:${PATH}
        unset PYTHONPATH
        unset PYTHONHOME
        ;;

    # for using conda installed Toil and Singularity to run CWL's
    toil)
        # need to get these env vars to propagate into the HPC jobs; SINGULARITY_DOCKER_USERNAME  SINGULARITY_DOCKER_PASSWORD
        [ -e /juno/work/ci/kellys5/projects/toil-settings/toil-settings.sh ] && \
        . /juno/work/ci/kellys5/projects/toil-settings/toil-settings.sh || \
        echo ">>> WARNING: could not find file toil-settings.sh, needed for SINGULARITY_DOCKER_USERNAME and SINGULARITY_DOCKER_PASSWORD; HPC jobs might break!"

        # NOTE: in the future might also need these
        # TOIL_LSF_ARGS
        # SINGULARITY_PULLDIR
        # SINGULARITY_CACHEDIR
        module load singularity/3.7.1
        export PATH=${PWD}/conda/bin:${PATH}
        export PATH="$(dirname $(which singularity))":${PATH}
        unset PYTHONPATH
        unset PYTHONHOME
        # need to unset this when running with Toil / cwltool ; need it to be set when doing 'singularity pull' commands
        unset SINGULARITY_PULLDIR
        ;;

    *)
        echo "unrecognized target called"
        # exit 1
        ;;
esac
