# pluto-cwl

**P**ost-processing & **L**ightweight **U**pdates **T**o pipeline **O**utput

CWL files and workflows to accompany the [helix_filters_01](https://github.com/mskcc/helix_filters_01) repo.

**NOTE:** see the help section under `make help` for the most up to date instructions.

## Test Suite

Development and testing takes place via the test suite.

Make sure that you have copies of the Singularity containers used cached locally

```
make singularity-pull
make singularity-pull-dev
make singularity-pull-facets
```


The included test suite can be run with:

```
make test
```

It typically takes about 45 minutes to run all included tests

- NOTE: tests require data sets that are pre-saved on the `juno` server

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

The run the script with the tests you are interested in, such as;

```
python tests/test_workflow_cwl.py
```

You can further select which test case(s) from the script you wish to run by adding their labels as args;

```
python tests/test_workflow_cwl.py TestClassName

python tests/test_workflow_cwl.py TestClassName.test_function
```

## Run a CWL

You can run a specific CWL workflow manually using the included `run.py` script. See `run.py -h` for available CWL's to run. 

Example usages;

- run a TMB workflow

```
$ ./run.py tmb_workflow --data-clinical examples/data_clinical.txt --assay-coverage 10000 --pairs examples/pairs.tsv
```

- run the main `workflow_with_facets` workflow

```
$ ./run.py workflow_with_facets --assay-coverage 100000 --project-id Project1 --cancer-type MEL --pairs examples/pairs.tsv --data-clinical examples/data_clinical.txt --sample-summary examples/sample_summary.txt --mutation-svs-txts examples/mutation_svs.txt --mutation-svs-mafs examples/mutation_svs_mafs.txt
```

NOTE: The `run.py` script still requires that the Singularity containers be present locally, see the section on running the test scripts for details. You might also have to set your environment using either `make bash` or `. env.juno.sh shell` first.
