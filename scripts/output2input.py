#!/usr/bin/env python3
"""
Script to parse the output of an Argos pipeline run
plus pairs file from output2pairs.py
and convert it into input format for pluto CLI operator

Usage
$ ./output2pairs.py /path/to/argos_output /path/to/argos_output/sample_pairing.txt pairs.tsv
$ ./output2input.py /path/to/argos_output /path/to/sample_data_clinical.txt /juno/work/ci/helix_filters_01/fixtures/demo/qc/demo_SampleSummary.txt pairs.tsv

-----
# parse the output of an Argos run
# need to gather up the tumor/normal files and output a formatted .tsv

assay_coverage
{
    "IMPACT341": 896637,
    "IMPACT468": 1139294,
    "IMPACT410": 1016335,
    "IMPACT505": 1213770
}


  --assay-coverage ASSAY_COVERAGE
  --project-id PROJECT_ID
  --cancer-type CANCER_TYPE
  --pairs PAIRS_FILE
  --data-clinical DATA_CLINICAL_FILE
  --sample-summary SAMPLE_SUMMARY_FILE
  --mutation-svs-txts MUTATION_SVS_TXT_FILES
  --mutation-svs-mafs MUTATION_SVS_MAF_FILES
  --cancer-study-identifier CANCER_STUDY_IDENTIFIER
  --project-name PROJECT_NAME
  --project-short-name PROJECT_SHORT_NAME
  --project-description PROJECT_DESCRIPTION
  --project-pi PROJECT_PI
  --request-pi REQUEST_PI
  --is-impact IS_IMPACT
  --argos-version-string ARGOS_VERSION_STRING
  --helix-filter-version HELIX_FILTER_VERSION
  --IMPACT-gene-list IMPACT_GENE_LIST
  --targets-list TARGETS_LIST
  --known-fusions KNOWN_FUSIONS_FILE
  --print-input
"""
import sys
import os
from collections import OrderedDict
import glob
import json

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import TableReader
from pluto.settings import DATA_SETS, KNOWN_FUSIONS_FILE, IMPACT_FILE
from operators.workflow_with_facets import WorkflowWithFacets, version
sys.path.pop(0)

# breadth of coverage for different target exome assays, used in TMB analysis
assay_coverages = {
    "IMPACT341": 896637,
    "IMPACT468": 1139294,
    "IMPACT410": 1016335,
    "IMPACT505": 1213770
}


def main():
    """
    """
    args = sys.argv[1:]
    argos_output_dir = args[0]
    sample_data_clinical_file = args[1] # sample_data_clinical_file = os.path.join(argos_output_dir, "sample_data_clinical.txt")
    sample_summary_file = args[2] # /juno/work/ci/helix_filters_01/fixtures/demo/qc/demo_SampleSummary.txt
    pairs_file = args[3] # the pairs.tsv that was output from the output2pairs.py script


    # default input dir locations
    bam_dir = os.path.join(argos_output_dir, "bam")
    maf_dir = os.path.join(argos_output_dir, "maf")

    # TODO: user input for these items
    project_id = "PROJECT_ID"
    cancer_type = "CANCER_TYPE"
    microsatellites_file = DATA_SETS['demo']['microsatellites_file']
    IMPACT_gene_list = IMPACT_FILE
    known_fusions_file = KNOWN_FUSIONS_FILE
    targets_list = DATA_SETS['demo']["targets_list"],

    # get the records from the data clinical file
    reader = TableReader(sample_data_clinical_file)
    data_clinical_records = [ rec for rec in reader.read() ]

    # get the pairs from the pairs file
    reader = TableReader(pairs_file)
    pairs_records = [ rec for rec in reader.read() ]

    # get the first assay from the list
    assay = [ rec["GENE_PANEL"] for rec in data_clinical_records ][0]
    try:
        assay_coverage = assay_coverages[assay]
    except KeyError: # try to parse the assay to get a value that matches
        assay = assay.split('+')[0] # IMPACT468+08390_Hg19
        assay_coverage = assay_coverages[assay]

    # get the lists of input maf files
    mutation_svs_txt_files = glob.glob(os.path.join(maf_dir, '*.svs.pass.vep.portal.txt'))
    mutation_svs_maf_files = glob.glob(os.path.join(maf_dir, '*.svs.pass.vep.maf'))
    normal_bam_files = []
    tumor_bam_files = []

    # get the bam files files
    for pair in pairs_records:
        # tumor_bam = glob.glob(os.path.join(bam_dir, "{}.*bam".format(pair["tumor_id"])))[0]
        # normal_bam = glob.glob(os.path.join(bam_dir, "{}.*bam".format(pair["normal_id"])))[0]
        normal_bam_files.append(pair["normal_bam"])
        tumor_bam_files.append(pair["tumor_bam"])

    operator = WorkflowWithFacets(
        pairs_file = pairs_file,
        assay_coverage = assay_coverage,
        project_id = project_id,
        cancer_type = cancer_type,
        data_clinical_file = sample_data_clinical_file,
        sample_summary_file = sample_summary_file,
        mutation_svs_txt_files = mutation_svs_txt_files,
        mutation_svs_maf_files = mutation_svs_maf_files,
        IMPACT_gene_list = IMPACT_gene_list,
        targets_list = targets_list,
        known_fusions_file = known_fusions_file,
        version = version,
        array_File_keys = None, # turn off reading in list of filepaths from .txt file for mutation_svs_txt_files, mutation_svs_maf_files
        list_File_keys = ['mutation_svs_txt_files', 'mutation_svs_maf_files'] # get file paths from lists instead of ^^
        )

    print(json.dumps(operator.input, indent = 4))


if __name__ == '__main__':
    main()
