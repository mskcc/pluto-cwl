# pluto-cwl

**P**ost-processing & **L**ightweight **U**pdates **T**o pipeline **O**utput

CWL files and workflows to accompany the [helix_filters_01](https://github.com/mskcc/helix_filters_01) repo.

**NOTE:** see the help section under `make help` for the most up to date instructions.

# Installation

To "install" the repo for running locally, use the following command;

```
make install singularity-pull-all
```

This will checkout the included `git` submodules, install a local `conda` with extra dependencies, and pull copies of the required Singularity containers locally.

To run manually in your current environment, source the environment script with the appropriate recipe;

```
# to run unit tests or run.py with cwltool
. env.juno.sh test
```

Alternatively, use the Makefile recipes since they already include environment configuration when running on Juno / Silo HPC servers.

# Usage

The primary entry point for the workflow is [`cwl/workflow_with_facets.cwl`](https://github.com/mskcc/pluto-cwl/blob/master/cwl/workflow_with_facets.cwl). 

## Run a CWL

### Setup

To run a CWL, you first need to make sure dependencies are available.

Running on Juno requires that the Singularity containers mentioned above are available;

```
make singularity-pull-all
```

If you are on Juno server using `cwltool` you can just load the included config for the `test` environment;

```
. env.juno.sh test
```

If you want to use Toil to submit jobs to LSF then you will need the `toil` environment config;

```
. env.juno.sh toil
```

### Run with `run.py`

-----

NOTE: using `run.py` is deprecated; instead check out [`run-cwltool.sh`](https://github.com/mskcc/pluto/blob/master/run-cwltool.sh) and [`run-toil.sh`](https://github.com/mskcc/pluto/blob/master/run-toil.sh)

-----

The included `run.py` script presents an easy command line interface to run CWL files, by mapping standard CLI flags and arguments to the CWL input arguments via the included [`operators`](https://github.com/mskcc/pluto-cwl/tree/master/operators). The script will automatically create the required `input.json` file for you, set up temporary directories for pipeline execution, and run the desired CWL with either the `cwltool` reference runner or Toil. See `run.py -h` for available options and CWL's to run.

Example usages;

- run a TMB workflow

```
$ ./run.py tmb_workflow --data-clinical examples/data_clinical.txt --assay-coverage 10000 --pairs examples/pairs.tsv
```

- run the main `workflow_with_facets` workflow

```
$ ./run.py workflow_with_facets --assay-coverage 100000 --project-id Project1 --cancer-type MEL --pairs examples/pairs.tsv --data-clinical examples/data_clinical.txt --sample-summary examples/sample_summary.txt --mutation-svs-txts examples/mutation_svs.txt --mutation-svs-mafs examples/mutation_svs_mafs.txt
```

Output will be in the `cwltool_output` or `toil_output` directories. Note that this includes `tmp` and `work` directories for the run, which may need to be deleted periodically.

### Run with `run-json.py`

If you already have a pre-made JSON file (such as debugging a production pipeline run) and want to run it easily, you can use the `run-json.py` script instead. 

## Test Suite

Development and testing takes place via the test suite.

Make sure that you have copies of the Singularity containers used cached locally

```
make singularity-pull-all
```

The included test suite can be run with:

```
make test
```

It typically takes about 45 minutes to run all included tests

- NOTE: tests require data sets that are pre-saved on the `juno` server

Some very large integration tests are skipped by default. To include all tests, export the environment variable `LARGE_TESTS=True` or include it in the command line invocation.

### Parallel Test Suite

An extra recipe is included which can run the tests in parallel, for example to run 8 tests at once you can use this command:

```
make test3 -j 8
```

### Single Test

For development purposes, it is helpful to be able to run only a specific test case, or subset of tests.

To do this, first enter an interactive bash session with the environment populated;

```
make bash
```

Or source the environment config file

```
. env.juno.sh test
```

Then run the script with the tests you are interested in, such as;

```
python tests/test_workflow_cwl.py
```

You can further select which test case(s) from the script you wish to run by adding their labels as args;

```
python tests/test_workflow_cwl.py TestClassName

python tests/test_workflow_cwl.py TestClassName.test_function
```


