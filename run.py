#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for running the CWL workflows
"""
import argparse
# from operators import workflow_with_facets
from operators.tmb_workflow import TMBWorkflow
from operators.workflow_with_facets import WorkflowWithFacets
from operators.example_workflow import ExampleWorkflow
from operators.env import EnvCWL
from operators.env_container import EnvContainerCWL
from operators.input import generate_sample_summary, generate_pairs_sheet, generate_data_clinical

def main():
    """
    Main function for CLI parsing
    """
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument("--engine", default = 'cwltool', dest = 'engine', choices = ['cwltool', 'toil'], help = "CWL execution engine to use")
    parser.add_argument("--print-command", action = 'store_true', dest = 'print_command', help = "Print the CWL runner command and exit")

    subparsers = parser.add_subparsers(help ='Sub-commands available', required = True)

    # Generate a blank template file
    _generate_data_clinical = subparsers.add_parser('data_clinical_template', help = 'Generate a blank template data clinical file in cBioPortal format')
    _generate_data_clinical.add_argument('--output', dest = 'output_file', default = "data_clinical.txt")
    _generate_data_clinical.set_defaults(func = generate_data_clinical)

    # Generate a blank template file
    _generate_pairs = subparsers.add_parser('pairs_template', help = 'Generate a blank template pairs file')
    _generate_pairs.add_argument('--output', dest = 'output_file', default = "pairs.tsv")
    _generate_pairs.set_defaults(func = generate_pairs_sheet)

    # Generate a blank template file
    _generate_sample_summary = subparsers.add_parser('sample_summary', help = 'Generate a blank template sample sumamry file')
    _generate_sample_summary.add_argument('--output', dest = 'output_file', default = "sample_summary.txt")
    _generate_sample_summary.set_defaults(func = generate_sample_summary)

    # Example workflow
    _example_workflow = subparsers.add_parser('example_workflow', help = 'Run the example workflow')
    _example_workflow.add_argument('--value', dest = 'value', required = True)
    _example_workflow.add_argument('--sampleIDs', dest = 'sampleIDs', nargs='+', required = True)
    _example_workflow.set_defaults(func = ExampleWorkflow._run)
    """
    $ ./run.py example_workflow --value foo --sampleIDs 1 2 3
    """

    # CWL runs for checking execution environment for debugging
    _env_cwl = subparsers.add_parser('env', help = 'Check the execution environment')
    _env_cwl.set_defaults(func = EnvCWL._run)

    _env_container_cwl = subparsers.add_parser('env_container', help = 'Check the execution environment inside the container')
    _env_container_cwl.set_defaults(func = EnvContainerCWL._run)

    # TMB workflow
    _tmb_workflow = subparsers.add_parser('tmb_workflow', help = 'Run the TMB workflow')
    _tmb_workflow.add_argument('--data-clinical', dest = 'data_clinical_file', required = True)
    _tmb_workflow.add_argument('--assay-coverage', dest = 'assay_coverage', required = True)
    _tmb_workflow.add_argument('--pairs', dest = 'pairs_file', required = True)
    _tmb_workflow.add_argument('--dir', dest = 'dir', help = 'Directory for pipeline execution and output')
    _tmb_workflow.set_defaults(func = TMBWorkflow._run)
    """
    $ ./run.py tmb_workflow --data-clinical examples/data_clinical.txt --assay-coverage 10000 --pairs examples/pairs.tsv
    """

    # Main full workflow with Facets
    _workflow_with_facets = subparsers.add_parser('workflow_with_facets', help = 'Run the full workflow with Facets')
    # required args
    _workflow_with_facets.add_argument('--assay-coverage', dest = 'assay_coverage', required = True)
    _workflow_with_facets.add_argument('--project-id', dest = 'project_id', required = True)
    _workflow_with_facets.add_argument('--cancer-type', dest = 'cancer_type', required = True)
    # required file samplesheets
    _workflow_with_facets.add_argument('--pairs', dest = 'pairs_file', required = True)
    _workflow_with_facets.add_argument('--data-clinical', dest = 'data_clinical_file', required = True)
    _workflow_with_facets.add_argument('--sample-summary', dest = 'sample_summary_file', required = True)
    # optional file samplesheets
    _workflow_with_facets.add_argument('--mutation-svs-txts', dest = 'mutation_svs_txt_files') # mutation_svs.txt; Sample1.Sample2.svs.pass.vep.portal.txt
    _workflow_with_facets.add_argument('--mutation-svs-mafs', dest = 'mutation_svs_maf_files') # mutation_svs_mafs.txt; Sample1.Sample2.svs.pass.vep.maf
    # optional args
    _workflow_with_facets.add_argument('--cancer-study-identifier', dest = 'cancer_study_identifier')
    _workflow_with_facets.add_argument('--project-name', dest = 'project_name')
    _workflow_with_facets.add_argument('--project-short-name', dest = 'project_short_name')
    _workflow_with_facets.add_argument('--project-description', dest = 'project_description')
    _workflow_with_facets.add_argument('--project-pi', dest = 'project_pi')
    _workflow_with_facets.add_argument('--request-pi', dest = 'request_pi')
    _workflow_with_facets.add_argument('--is-impact', dest = 'is_impact', default = True)
    _workflow_with_facets.add_argument('--argos-version-string', dest = 'argos_version_string', default = '2.x')
    _workflow_with_facets.add_argument('--helix-filter-version', dest = 'helix_filter_version')
    _workflow_with_facets.add_argument('--IMPACT-gene-list', dest = 'IMPACT_gene_list')
    _workflow_with_facets.add_argument('--targets-list', dest = 'targets_list')
    _workflow_with_facets.add_argument('--known-fusions', dest = 'known_fusions_file')
    _workflow_with_facets.add_argument('--print-input', dest = 'print_input', action = 'store_true')
    _workflow_with_facets.set_defaults(func = WorkflowWithFacets._run)
    """
    $ ./run.py workflow_with_facets --assay-coverage 100000 --project-id Project1 --cancer-type MEL --pairs examples/pairs.tsv --data-clinical examples/data_clinical.txt --sample-summary examples/sample_summary.txt --mutation-svs-txts examples/mutation_svs.txt --mutation-svs-mafs examples/mutation_svs_mafs.txt --print-input
    """

    args = parser.parse_args()
    args.func(**vars(args))

if __name__ == '__main__':
    main()
