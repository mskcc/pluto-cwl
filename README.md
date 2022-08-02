# pluto-cwl

**P**ost-processing & **L**ightweight **U**pdates **T**o pipeline **O**utput

CWL files and workflows to accompany the [helix_filters_01](https://github.com/mskcc/helix_filters_01) repo. Supported by infrastructure in the [`pluto`](https://github.com/mskcc/pluto) submodule.

# Installation & Setup

Install dependencies for the repo with the command:

```
make install 
```

This will checkout the included `git` submodules and install a local `conda` with extra dependencies.

Use this command to activate the installed environment for running workflows:

```
. env.juno.sh toil
```

This will:

- update your environment to use the `cwltool` and `toil` installed in the local `conda`
- (if running on Juno HPC) update your environment with Toil variables needed to run on Juno
- (if running on Juno HPC) upate your environment to use pre-cached Singularity containers located on Juno


# Run a CWL

The primary entry point for the workflow is [`cwl/workflow_with_facets.cwl`](https://github.com/mskcc/pluto-cwl/blob/master/cwl/workflow_with_facets.cwl). 

You can run a CWL included in this repo by using the wrapper scripts bundled in the `pluto` submodule; 

- [`pluto/run-cwltool.sh`](https://github.com/mskcc/pluto/blob/master/run-cwltool.sh) for simple use cases
- [`pluto/run-toil.sh`](https://github.com/mskcc/pluto/blob/master/run-toil.sh) if parallel processing and HPC (LSF) useage is required

## Test Suite

Development and testing takes place via the test suite.

The included test suite can be run with:

```
make test
```

It typically takes about 45 minutes to run all included tests

- NOTE: tests require data sets that are pre-saved on the `juno` server

Some very large integration tests are skipped by default. To include all tests, export the environment variable `LARGE_TESTS=True` or include it in the command line invocation. You can also change the CWL engine from `cwltool` to `toil`, among other settings, the same way. For example;

```
LARGE_TESTS=True CWL_ENGINE=Toil PRINT_COMMAND=True TMP_DIR=/scratch USE_LSF=True make test
```

Available environment variable settings are derived from the [`pluto.settings`](https://github.com/mskcc/pluto/blob/master/settings.py) submodule.

### Parallel Test Suite

An extra recipe is included which can run the tests in parallel, for example to run 8 tests at once you can use this command:

```
make parallel-test
```

### Single Test

For development purposes, it is helpful to be able to run only a specific test case, or subset of tests.

You can run just the script with the tests you are interested in, such as;

```
python tests/test_workflow_cwl.py
```

You can further select which test case(s) from the script you wish to run by adding their labels as args;

```
python tests/test_workflow_cwl.py TestClassName

python tests/test_workflow_cwl.py TestClassName.test_function
```

This can be combined with the environment variables described above (such as `LARGE_TESTS`, `PRINT_COMMAND`, `KEEP_TMP`, etc.).
 
