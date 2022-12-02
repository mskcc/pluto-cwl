export SHELL:=/bin/bash
.ONESHELL:
# enforce some stricter shell options to
OLDSHELLOPTS:=$(SHELLOPTS)
export SHELLOPTS:=$(if $(SHELLOPTS),$(SHELLOPTS):)pipefail:errexit
UNAME:=$(shell uname)

# .sh file to source for each recipe
ENVSH:=env.juno.sh

define help
This is the Makefile for helix filters CWL

Development Guide
-----------------

Create a new .cwl file in the cwl/ dir, and a corresponding test script in the tests/ dir

Pull copies of the Singularity containers locally:

```
# pulls the latest tagged container
make singularity-pull

# pulls the latest master branch container
make singularity-pull-dev
```

Run the individual test case & cwl you are developing directly

```
# load the environment for running pluto-cwl
make bash

# run your test script
python3 tests/test_tmb_cwl.py
```

Note that `unittest.main()` needs to be included in your test script to run it directly from CLI

When you are done working on your updates, run the full test suite with:

```
# serial test runner
make test

# or

# run tests in parallel
make test3 -j 8

# or

# run integration tests only
make integration_test
```

To do a version bump of all CWL files to a new version of the helix_filters container, use:

```
make update-container-tags

# example
make update-container-tags OLD_TAG='20.11.2' NEW_TAG='21.01.0'
```



NOTE: these may be deprecated;
Run the workflow against the demo dataset with:

```
make run
```

Run the Facets Suite workflow against the demo dataset with
```
make facets
```
endef
export help
help:
	@printf "$$help"
.PHONY : help



# ~~~~~ Install Dependencies ~~~~~~ #
# This is only needed if you are not running on Juno/Silo server, or want to use Toil
# NOTE: see env.sh for PATH and PYTHONHOME, PYTHONPATH modifications for this to work correctly
ifeq ($(UNAME), Darwin)
CONDASH:=Miniconda3-4.7.12.1-MacOSX-x86_64.sh
endif

ifeq ($(UNAME), Linux)
CONDASH:=Miniconda3-4.7.12.1-Linux-x86_64.sh
endif

CONDAURL:=https://repo.anaconda.com/miniconda/$(CONDASH)

conda: SHELLOPTS=$(OLDSHELLOPTS)
conda:
	@echo ">>> Setting up conda..."
	@wget "$(CONDAURL)" && \
	bash "$(CONDASH)" -b -p conda && \
	rm -f "$(CONDASH)"

# source conda/bin/activate
# conda deactivate
install: SHELLOPTS=$(OLDSHELLOPTS)
install: conda
	. "$(ENVSH)" conda && \
	. conda/bin/activate && \
	conda env update --file environment.yml

# conda install -y conda-forge::jq=1.5
# pip install -r requirements.txt
# $(MAKE) init

init: SHELLOPTS=$(OLDSHELLOPTS)
init:
	git submodule init
	git submodule update






# ~~~~~ Container ~~~~~ #
# pull the Docker container and convert it to Singularity container image file
export SINGULARITY_CACHEDIR:=/juno/work/ci/pluto-cwl-test/cache
# GIT_TAG:=$(shell git describe --tags --abbrev=0)
HF_CONTAINER:=mskcc/helix_filters_01
HF_TAG:=21.3.5
DOCKER_TAG:=$(HF_CONTAINER):$(HF_TAG)
DOCKER_DEV_TAG:=$(HF_CONTAINER):latest
# NOTE: you cannot use a filename with a ':' as a Makefile target
SINGULARITY_SIF:=mskcc_helix_filters_01:$(HF_TAG).sif
SINGULARITY_DEV_SIF:=mskcc_helix_filters_01:latest.sif
singularity-pull:
	. "$(ENVSH)" toil && \
	singularity pull --force --name "$(SINGULARITY_SIF)" docker://$(DOCKER_TAG)

singularity-pull-dev:
	. "$(ENVSH)" toil && \
	singularity pull --force --name "$(SINGULARITY_DEV_SIF)" docker://$(DOCKER_DEV_TAG)
# unset SINGULARITY_CACHEDIR && \
# module load singularity/3.3.0 && \

# shell into the Singularity container to check that it looks right
singularity-shell:
	- . "$(ENVSH)" toil && \
	singularity shell "$(SINGULARITY_SIF)"

# mskcc/roslin-variant-facets:1.6.3
# mskcc/helix_filters_01:facets-1.6.3
FACETS_DOCKERTAG:=mskcc/helix_filters_01:facets-1.6.3
FACETS_SIF:=mskcc_helix_filters_01:facets-1.6.3.sif
singularity-pull-facets:
	. "$(ENVSH)" toil && \
	singularity pull --force --name "$(FACETS_SIF)" docker://$(FACETS_DOCKERTAG)

FACETS_SUITE_DOCKERTAG:=mskcc/helix_filters_01:facets-suite-2.0.6
FACETS_SUITE_SIF:=helix_filters_01_facets-suite-2.0.6.sif
singularity-pull-facets-suite:
	. "$(ENVSH)" toil && \
	singularity pull --force --name "$(FACETS_SUITE_SIF)" docker://$(FACETS_SUITE_DOCKERTAG)

# docker://cmopipeline/getbasecountsmultisample:1.2.2
FILLOUT_DOCKERTAG:=mskcc/helix_filters_01:getbasecountsmultisample-1.2.2
FILLOUT_SIF:=helix_filters_01_getbasecountsmultisample-1.2.2.sif
singularity-pull-fillout:
	. "$(ENVSH)" toil && \
	singularity pull --force --name "$(FILLOUT_SIF)" docker://$(FILLOUT_DOCKERTAG)

MSI_DOCKERTAG:=mskcc/msisensor:0.2
MSI_SIF:=mskcc_msisensor:0.2.sif
singularity-pull-msi:
	. "$(ENVSH)" toil && \
	singularity pull --force --name "$(MSI_SIF)" docker://$(MSI_DOCKERTAG)

IGV_REPORTS_DOCKERTAG:=$(HF_CONTAINER):igv-reports-1.0.1
IGV_REPORTS_SIF:=mskcc_helix_filters_01:igv-reports-1.0.1.sif
singularity-pull-igv-reports:
	. "$(ENVSH)" toil && \
	singularity pull --force --name "$(IGV_REPORTS_SIF)" docker://$(IGV_REPORTS_DOCKERTAG)

# mskcc/roslin-variant-cmo-utils:1.9.15
# mskcc_roslin-variant-cmo-utils:1.9.15.sif
CMOUTILS_DOCKERTAG:=mskcc/roslin-variant-cmo-utils:1.9.15
CMOUTILS_SIF:=mskcc_roslin-variant-cmo-utils:1.9.15.sif
singularity-pull-cmoutils:
	. "$(ENVSH)" toil && \
	singularity pull --force --name "$(CMOUTILS_SIF)" docker://$(CMOUTILS_DOCKERTAG)

R_DOCKERTAG:=mskcc/helix_filters_01:R-3.5.1
R_SIF:=mskcc_helix_filters_01:R-3.5.1.sif
singularity-pull-r:
	. "$(ENVSH)" toil && \
	singularity pull --force --name "$(R_SIF)" docker://$(R_DOCKERTAG)

REPORT_DOCKERTAG:=mskcc/helix_filters_01:reporting
REPORT_SIF:=mskcc_helix_filters_01:reporting.sif
singularity-pull-report:
	. "$(ENVSH)" toil && \
	singularity pull --force --name "$(REPORT_SIF)" docker://$(REPORT_DOCKERTAG)


# Need this container for the VEP cache dir; its too large to pass as CWL workflow inputs
# https://hub.docker.com/r/mskcc/roslin-variant-vcf2maf/tags?page=1&ordering=last_updated
# if there is a local copy of the file then copy it here
VEP_SIF:=roslin-variant-vcf2maf_1.6.17.sif
VEP_SIF_LOCAL:=/juno/work/ci/singularity_images/roslin-variant-vcf2maf_1.6.17.sif
VEP_DOCKERTAG:=mskcc/roslin-variant-vcf2maf:1.6.17
singularity-pull-vep:
	if [ ! -e "$(VEP_SIF)" ]; then
	if [ -e "$(VEP_SIF_LOCAL)" ]; then
	ln -s "$(VEP_SIF_LOCAL)" "$(VEP_SIF)"
	else
	. "$(ENVSH)" toil && singularity pull --force --name "$(VEP_SIF)" docker://$(VEP_DOCKERTAG)
	fi
	fi
# rsync -vrthP "$(VEP_SIF_LOCAL)" "$(VEP_SIF)"

singularity-pull-all: singularity-pull singularity-pull-dev singularity-pull-facets-suite singularity-pull-facets singularity-pull-fillout singularity-pull-igv-reports singularity-pull-cmoutils singularity-pull-msi singularity-pull-vep singularity-pull-r

# change the Docker tag for all the CWL files from the old pattern to the new pattern
OLD_TAG:=
NEW_TAG:=
UPDATE_CONTAINER:=$(HF_CONTAINER)
update-container-tags:
	[ -z "$(OLD_TAG)" ] && echo "ERROR: OLD_TAG value missing" && exit 1 || :
	[ -z "$(NEW_TAG)" ] && echo "ERROR: NEW_TAG value missing" && exit 1 || :
	for i in $$(find cwl -type f -exec grep -l 'dockerPull: $(UPDATE_CONTAINER)' {} \;); do \
	perl -i -pe 's/$(OLD_TAG)/$(NEW_TAG)/g' $$i ; \
	done







# ~~~~~ Debug & Development ~~~~~ #
# run all individual test cases in parallel

# RECIPES:
# parallel-test; run all tests in parallel
# parallel-test-log; run all tests in parallel but save stdout to log file (stdout might not appear in console right away but its running! check htop if you arent sure if its running)


# EXAMPLES:
# initialize the "toil" environment (assuming you have installed dependencies)
# $ . env.juno.sh toil
# run all tests with default settings
# $ make parallel-test

# run with logging
# $make parallel-test-log

# run only test from one .py file, run only 4 tests in parallel
# $make parallel-test-log T=tests/test_add_af_cwl.py P=4

# pass env vars to control test execution
# $ PRINT_TESTNAME=T make parallel-test T=tests/test_generate_cBioPortal_file_cwl.py
# $ CWL_ENGINE=Toil TMP_DIR=/scratch PRINT_COMMAND=T PRINT_TESTNAME=True make parallel-test P=20
# $ CWL_ENGINE=Toil TMP_DIR=/fscratch/tmp PRINT_COMMAND=T PRINT_TESTNAME=T make parallel-test
# $ CWL_ENGINE=Toil TMP_DIR=/fscratch/tmp LARGE_TESTS=True PRINT_COMMAND=True PRINT_TESTNAME=True INTEGRATION_TESTS=True make parallel-test-log
# $ CWL_ENGINE=Toil USE_LSF=T PRINT_COMMAND=T PRINT_TESTNAME=T make parallel-test-log

# NOTE: ^^^ env vars must come BEFORE `make ....`, and Makefile vars must come AFTER

# number of parallel tasks;
P:=16
# test target; can also be an individual test_*.py file!
T:=tests/
# NOTE: the logging here can cause delay before lines are printed to file or console ; stdbuf -oL ; this breaks some R sigularity containers
_LOGDATE:=$(shell date +%s)
_LOGDATE_LONG:=$(shell date)
_LOGFILE:=testing.$(_LOGDATE).log
parallel-test-log:
	echo ">>>>> ------ start parallel-test-log $(_LOGDATE_LONG) ($(_LOGDATE)) -------- <<<<<" >> "$(_LOGFILE)"
	time { QUIET=True ./print_tests.py "$(T)" | \
	xargs -n 1 -P "$(P)" nice python3 -m unittest ; } 2>&1 | \
	tee -a "$(_LOGFILE)"
	echo ">>>>> ------ stop parallel-test-log $(_LOGDATE_LONG) ($(_LOGDATE)) -------- <<<<<" >> "$(_LOGFILE)"


# same but without logging
parallel-test:
	./print_tests.py "$(T)" | \
	xargs -n 1 -P "$(P)" nice python3 -m unittest


pytest:
	. "$(ENVSH)" toil && \
	source conda/bin/activate && \
	nice pytest -n auto --maxprocesses 24 tests



#
# OLD TEST RUNNING RECIPES:
#

# Run the test suite
# NOTE: run with `$ LARGE_TESTS=True python3 tests/... ` to enable large test cases

# TODO: why is this here and do we need it?? might be used in pluto settings.py??
export FIXTURES_DIR:=/juno/work/ci/helix_filters_01/fixtures

# TODO: figure out if we can run the tests in parallel or otherwise make it faster
# for some reason the test recipe is not running all tests....
test:
	. "$(ENVSH)" toil && \
	if [ ! -e "$(SINGULARITY_SIF)" ]; then $(MAKE) singularity-pull; fi && \
	for i in tests/test_*.py; do echo $$i; python3 $$i; done

# run tests in parallel;
# $ make test3 -j 4
# NOTE: this runs each test_*.py file in parallel; does NOT run individual test cases in parallel!!
TEST_ENV:=toil
TESTS:=$(shell ls tests/test_*.py)
$(TESTS):
	. "$(ENVSH)" "$(TEST_ENV)" && echo $@; python3 $@
.PHONY:$(TESTS)
test3:$(TESTS)


# test recipe for Jenkins
# TODO: review this, it has problems with the pwd during Jenkins execution; consider only running medium or XL test
integration_test:
	. "$(ENVSH)" integration_test && \
	cd pluto && \
	python test_tools.py && \
	python test_serializer.py && \
	cd .. && \
	#for i in tests/test_*workflow*.py; do echo $$i; python3 $$i; rm -rf $TMP_DIR/tmp* /scratch/jenkins/tmp*; done
	python tests/test_workflow_with_facets.xl.py && \
	python tests/test_workflow_with_facets.medium.py

# alternative Jenkins test recipe
# TODO: do we still need this?
integration_test_1:
	. "$(ENVSH)" integration_test && \
	for i in tests/test_*cwl*.py; do echo $$i; python3 $$i; done
