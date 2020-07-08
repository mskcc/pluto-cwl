# pluto-cwl

**P**ost-processing & **L**ightweight **U**pdates **T**o pipeline **O**utput

CWL files and workflows to accompany the *pluto* repo (in development)

## Usage

The main workflows can be run from the included `Makefile`. These recipes are set up for usage on MSKCC `juno` and `silo` servers and might need some modification for usage elsewhere. 

First, clone this repo;

```
git clone https://github.com/mskcc/pluto-cwl.git
cd pluto-cwl
```

The most recent instructions can be seen with

```
make help
```

The main workflow can be run against a test dataset saved on the `juno` server with the command 

```
make run
```

Extra arguments can be supplied to direct the workflow to run against other datasets (see `Makefile` for details), example;

```
make run \
PROJ_ID=My_Project \
TARGETS_LIST=/path/to/targest.ilist \
DATA_DIR=/path/to/pipeline_data/
```

Depending on your user account settings on `juno` you might need to pull a copy of the needed Singularity containers manually before running;

```
make singularity-pull
make singularity-pull-facets
```

## Run the Test Suite

The included test suite can be run with:

```
make test
```

- NOTE: tests require data sets that are pre-saved on the `juno` server

## Related

`pluto-cwl` is a migration of the CWL workflows that were originally included with the `helix_filters_01` repository; https://github.com/mskcc/helix_filters_01
