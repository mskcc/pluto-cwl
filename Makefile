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

install: SHELLOPTS=$(OLDSHELLOPTS)
install: conda
	. "$(ENVSH)" conda
	conda install -y conda-forge::jq=1.5
	pip install -r requirements.txt

# ~~~~~ Container ~~~~~ #
# pull the Docker container and convert it to Singularity container image file
export SINGULARITY_CACHEDIR:=/juno/work/ci/pluto-cwl-test/cache
# GIT_TAG:=$(shell git describe --tags --abbrev=0)
DOCKER_TAG:=mskcc/helix_filters_01:21.01.1
DOCKER_DEV_TAG:=mskcc/helix_filters_01:latest
# NOTE: you cannot use a filename with a ':' as a Makefile target
SINGULARITY_SIF:=mskcc_helix_filters_01:21.01.1.sif
SINGULARITY_DEV_SIF:=mskcc_helix_filters_01:latest.sif
singularity-pull:
	. "$(ENVSH)" singularity && \
	singularity pull --force --name "$(SINGULARITY_SIF)" docker://$(DOCKER_TAG)

singularity-pull-dev:
	. "$(ENVSH)" singularity && \
	singularity pull --force --name "$(SINGULARITY_DEV_SIF)" docker://$(DOCKER_DEV_TAG)
# unset SINGULARITY_CACHEDIR && \
# module load singularity/3.3.0 && \

# shell into the Singularity container to check that it looks right
singularity-shell:
	- . "$(ENVSH)" singularity && \
	singularity shell "$(SINGULARITY_SIF)"

FACETS_DOCKERTAG:=stevekm/facets-suite:2.0.6
FACETS_SIF:=stevekm_facets-suite:2.0.6.sif
singularity-pull-facets:
	. "$(ENVSH)" singularity && \
	singularity pull --force --name "$(FACETS_SIF)" docker://$(FACETS_DOCKERTAG)

# change the Docker tag for all the CWL files from the old pattern to the new pattern
OLD_TAG:=20.06.1
NEW_TAG:=20.06.2
update-container-tags:
	for i in $$(find cwl -type f -exec grep -l 'dockerPull: mskcc/helix_filters_01' {} \;); do \
	perl -i -pe 's/$(OLD_TAG)/$(NEW_TAG)/g' $$i ; \
	done


# ~~~~~ Debug & Development ~~~~~ #
# Run the test suite
export FIXTURES_DIR:=/juno/work/ci/helix_filters_01/fixtures
# TODO: figure out why this is missing some tests
test2:
	. "$(ENVSH)" test && \
	if [ ! -e "$(SINGULARITY_SIF)" ]; then $(MAKE) singularity-pull; fi && \
	python3 test.py

# TODO: figure out if we can run the tests in parallel or otherwise make it faster
# for some reason the test recipe is not running all tests....
test:
	. "$(ENVSH)" test && \
	if [ ! -e "$(SINGULARITY_SIF)" ]; then $(MAKE) singularity-pull; fi && \
	for i in tests/test_*.py; do echo $$i; python3 $$i; done

# run tests in parallel;
# $ make test3 -j 4
TESTS:=$(shell ls tests/test_*.py)
$(TESTS):
	. "$(ENVSH)" test && echo $@; python3 $@
.PHONY:$(TESTS)
test3:$(TESTS)


# interactive session with environment populated
bash: ENV=shell
bash:
	. "$(ENVSH)" "$(ENV)" && bash

clean:
	rm -rf cache tmp

clean-all: clean
	rm -rf output portal analysis mutation_maf_files.txt facets_hisens_seg_files.txt facets_hisens_cncf_files.txt mutation_svs_txt_files.txt mutation_svs_maf_files.txt


# Example args for running with Toil
run-toil: OUTPUTDIR=$(CURDIR)/toil_output
run-toil: LOGFILE=$(OUTPUTDIR)/toil.log
run-toil: JOBSTORE=$(OUTPUTDIR)/job-store
run-toil: WORKDIR=$(OUTPUTDIR)/work
run-toil: ENV=toil
run-toil:
	. $(ENVSH) $(ENV)
	unset SINGULARITY_CACHEDIR
	mkdir -p "$(OUTPUTDIR)"
	mkdir -p "$(WORKDIR)"
	[ -e "$(JOBSTORE)" ] && rm -rf "$(JOBSTORE)" || :
	( toil-cwl-runner \
	--logFile "$(LOGFILE)" \
	--outdir "$(OUTPUTDIR)" \
	--workDir "$(WORKDIR)" \
	--jobStore "$(JOBSTORE)" \
	--singularity \
	--batchSystem lsf --disableCaching True \
	--disable-user-provenance --disable-host-provenance \
	--preserve-entire-environment \
	cwl/example_workflow.cwl cwl/example_input.json ) > toil.stdout.txt
# --preserve-environment PATH TMPDIR TOIL_LSF_ARGS SINGULARITY_CACHEDIR SINGULARITY_TMPDIR SINGULARITY_PULLDIR PWD \
# --maxLocalJobs 500 \
