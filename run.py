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
from operators.ls import LsCWL
from operators.ls_dir import LsDirCWL
from operators.concat_tables_dir import ConcatTablesDirCWL
from operators.consensus_bed import ConsensusBed
from operators.concat_mafs import ConcatMafs
from operators.deduplicate_maf import DeduplicateMaf
from operators.consensus_maf import ConsensusMaf
from operators.samples_fillout import SamplesFillout
from operators.maf_filter import MafFilter
from operators.maf2vcf import Maf2Vcf
from operators.vcf_sort import VcfSort
from operators.bgzip import Bgzip
from operators.index_vcf import IndexVcf
from operators.igv_reports import IgvReports
from operators.run_facets import RunFacets
from operators.snp_pileup import SnpPileup
from operators.head import HeadCWL
from operators.input import generate_sample_summary, generate_pairs_sheet, generate_data_clinical, generate_samples_fillout_sheet

def main():
    """
    Main function for CLI parsing
    """
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument("--engine", default = 'cwltool', dest = 'engine', choices = ['cwltool', 'toil'], help = "CWL execution engine to use")
    parser.add_argument("--print-command", action = 'store_true', dest = 'print_command', help = "Print the CWL runner command and exit")
    parser.add_argument("--restart", action = 'store_true', dest = 'restart', help = "Restart a previous run; requires jobStore")
    parser.add_argument("--debug", action = 'store_true', dest = 'debug', help = "Restart a previous run; requires jobStore")
    parser.add_argument("--jobStore", dest = 'jobStore', default = None, help = "Job store to use for a restarted run")
    # parser.add_argument("--cleanWorkDir", dest = 'cleanWorkDir', default = 'onSuccess', help = "When to clean the work dir")
    # parser.add_argument("--clean", dest = 'clean', default = 'onSuccess', help = "When to clean the work dir")

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

    _generate_samples_fillout_sheet = subparsers.add_parser('samples_fillout_sheet', help = 'Generate a blank template sample fillout samplesheet file')
    _generate_samples_fillout_sheet.add_argument('--output', dest = 'output_file', default = "samples.fillout.tsv")
    _generate_samples_fillout_sheet.set_defaults(func = generate_samples_fillout_sheet)



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

    _ls_cwl = subparsers.add_parser('ls', help = 'Check the execution directory')
    _ls_cwl.set_defaults(func = LsCWL._run)

    _ls_dir_cwl = subparsers.add_parser('ls_dir', help = 'Check the execution directory with files input')
    _ls_dir_cwl.add_argument('input_files', nargs='+', help="Input files to stage in the directory")
    _ls_dir_cwl.set_defaults(func = LsDirCWL._run)

    _head_cwl = subparsers.add_parser('head', help = 'Head a file')
    _head_cwl.add_argument('input_file', help="Input files to stage in the directory")
    _head_cwl.add_argument('--num-lines', dest = 'num_lines', required = True)
    _head_cwl.set_defaults(func = HeadCWL._run)



    _concat_tables_dir = subparsers.add_parser('concat_tables_dir', help = 'Concatenate tables')
    _concat_tables_dir.add_argument('input_files', nargs='*', help="Input files to stage in the directory")
    _concat_tables_dir.add_argument('--input-files-list', dest = 'input_files_list', help="List of input files")
    _concat_tables_dir.add_argument('--comments', dest = 'comments', action = "store_true", help="Parse file header comments")
    _concat_tables_dir.add_argument('--na-str', dest = 'na_str', default = ".", help="NA string")
    _concat_tables_dir.add_argument('-o', '--output-file', dest = 'output_filename', default = "output.txt", help="Output filename")
    _concat_tables_dir.set_defaults(func = ConcatTablesDirCWL._run)

    _concat_tables_dir = subparsers.add_parser('concat_mafs', help = 'Concatenate maf files')
    _concat_tables_dir.add_argument('input_files', nargs='*', help="Input maf files")
    _concat_tables_dir.add_argument('--input-files-list', dest = 'input_files_list', help="List of input files")
    _concat_tables_dir.add_argument('-o', '--output-file', dest = 'output_filename', default = "output.maf", help="Output filename")
    _concat_tables_dir.set_defaults(func = ConcatMafs._run)

    _concat_tables_dir = subparsers.add_parser('deduplicate_maf', help = 'Deduplicate and sort rows in a maf file')
    _concat_tables_dir.add_argument('input_file', help="Input maf files")
    # _concat_tables_dir.add_argument('--input-files-list', dest = 'input_files_list', help="List of input files")
    # _concat_tables_dir.add_argument('-o', '--output-file', dest = 'output_filename', default = "output.maf", help="Output filename")
    _concat_tables_dir.set_defaults(func = DeduplicateMaf._run)

    _consensus_bed = subparsers.add_parser('consensus_bed', help = 'Merge maf files into a consensus bed file')
    _consensus_bed.add_argument('maf_files', nargs='*', help="Input maf files")
    _consensus_bed.add_argument('--maf-files-list', dest = 'maf_files_list', help="List of input files")
    _consensus_bed.set_defaults(func = ConsensusBed._run)

    _consensus_maf = subparsers.add_parser('consensus_maf', help = 'Merge maf files into a consensus maf file')
    _consensus_maf.add_argument('maf_files', nargs='*', help="Input maf files")
    _consensus_maf.add_argument('--maf-files-list', dest = 'maf_files_list', help="List of input files")
    _consensus_maf.set_defaults(func = ConsensusMaf._run)

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

    _samples_fillout = subparsers.add_parser('samples_fillout', help = 'Run the samples fillout workflow')
    _samples_fillout.add_argument('--samplesheet', dest = 'samplesheet', required = True)
    _samples_fillout.add_argument('--ref-fasta', dest = 'ref_fasta', required = True)
    _samples_fillout.set_defaults(func = SamplesFillout._run)
    """
    $ ./run.py samples_fillout --samplesheet examples/samples.fillout.tsv --ref-fasta /juno/work/ci/resources/genomes/GRCh37/fasta/b37.fasta
    """

    _maf_filter = subparsers.add_parser('maf_filter', help = 'Run the maf filter workflow')
    _maf_filter.add_argument('--maf-file', dest = 'maf_file', required = True)
    _maf_filter.add_argument('--argos-version-string', dest = 'argos_version_string', default = "2.x")
    _maf_filter.add_argument('--is-impact', dest = 'is_impact', action = "store_false")
    _maf_filter.add_argument('--analysis-mutations-filename', dest = 'analysis_mutations_filename', default = "analysis.muts.maf")
    _maf_filter.add_argument('--cbio-mutation-data-filename', dest = 'cbio_mutation_data_filename', default = "data_mutations_extended.txt")
    _maf_filter.add_argument('--rejected-file', dest = 'rejected_file', default = "rejected.muts.maf")
    _maf_filter.add_argument('--keep-rejects', dest = 'keep_rejects', action = "store_false")
    _maf_filter.set_defaults(func = MafFilter._run)
    """
    $ ./run.py maf_filter --maf-file /juno/work/ci/helix_filters_01/fixtures/Proj_08390_G/maf/Sample6.Sample5.muts.maf
    """


    _maf2vcf = subparsers.add_parser('maf2vcf', help = 'Run the maf2vcf workflow')
    _maf2vcf.add_argument('--maf-file', dest = 'maf_file', required = True)
    _maf2vcf.add_argument('--ref-fasta', dest = 'ref_fasta', required = True)
    _maf2vcf.add_argument('--output-filename', dest = 'output_vcf_filename', default = "output.vcf")
    _maf2vcf.set_defaults(func = Maf2Vcf._run)
    """
    $ ./run.py maf2vcf --maf-file cwltool_output//output/analysis.muts.maf --ref-fasta /juno/work/ci/resources/genomes/GRCh37/fasta/b37.fasta
    """

    _vcf_sort = subparsers.add_parser('vcf_sort', help = 'Run the vcf_sort workflow')
    _vcf_sort.add_argument('--vcf-file', dest = 'vcf_file', required = True)
    _vcf_sort.add_argument('--output-filename', dest = 'output_filename', default = "output.sorted.vcf")
    _vcf_sort.set_defaults(func = VcfSort._run)
    """
    $ ./run.py vcf_sort --vcf-file cwltool_output/output/output.vcf
    """

    _bgzip = subparsers.add_parser('bgzip', help = 'Run the bgzip workflow')
    _bgzip.add_argument(dest = 'input_file')
    _bgzip.add_argument('--output-filename', dest = 'output_filename', default = "output.gz")
    _bgzip.set_defaults(func = Bgzip._run)
    """
    $ ./run.py bgzip cwltool_output/output/output.sorted.vcf --output-filename output.sorted.vcf.gz
    """

    _index_vcf = subparsers.add_parser('index_vcf', help = 'Run the index_vcf workflow')
    _index_vcf.add_argument(dest = 'input_file')
    _index_vcf.set_defaults(func = IndexVcf._run)
    """
    $ ./run.py index_vcf cwltool_output/output/output.sorted.vcf.gz
    """

    _igv_report = subparsers.add_parser('igv_report', help = 'Run the index_vcf workflow')
    _igv_report.add_argument('--ref-fasta', dest = 'ref_fasta', required = True)
    _igv_report.add_argument('--sites', dest = 'sites', required = True)
    _igv_report.add_argument('--vcf-files', dest = 'vcf_gz_files', nargs='+', help="Input vcf.gz files. Should have adjacent .tbi files")
    _igv_report.add_argument('--bam-files', dest = 'bam_files', nargs='+', help="Input bam files. Should have adjacent .bai files")
    _igv_report.add_argument('--output-filename', dest = 'output_filename', default = "igv.html")
    _igv_report.set_defaults(func = IgvReports._run)
    """
    $ ./run.py igv_report --ref-fasta /juno/work/ci/resources/genomes/GRCh37/fasta/b37.fasta --vcf-files cwltool_output/output/output.sorted.vcf.gz --bam-files /juno/work/ci/helix_filters_01/fixtures/Proj_08390_G/bam/Sample5.rg.md.abra.printreads.bam /juno/work/ci/helix_filters_01/fixtures/Proj_08390_G/bam/Sample6.rg.md.abra.printreads.bam --sites cwltool_output/output/output.sorted.vcf.gz
    """

    _run_facets = subparsers.add_parser('run_facets', help = 'Run the Facets wrapper workflow')
    _run_facets.add_argument('--snp-pileup', dest = 'snp_pileup', required = True)
    _run_facets.add_argument('--sample-id', dest = 'sample_id', required = True)
    _run_facets.add_argument('--purity-cval', dest = 'purity_cval', default = "100")
    _run_facets.add_argument('--cval', dest = 'cval', default = "50")
    _run_facets.add_argument('--seed', dest = 'seed', default = "1000")
    _run_facets.add_argument('--min-nhet', dest = 'min_nhet', default = "25")
    _run_facets.add_argument('--purity-min-nhet', dest = 'purity_min_nhet', default = "25")
    _run_facets.set_defaults(func = RunFacets._run)
    """
    $ ./run.py run_facets --snp-pileup /juno/work/ci/helix_filters_01/fixtures/Proj_08390_G/snp-pileup/Sample34.Sample33.snp_pileup.gz --sample-id Sample34.Sample33
    """

    _snp_pileup = subparsers.add_parser('snp_pileup', help = 'Run the Facets wrapper workflow')
    _snp_pileup.add_argument('--snps-vcf', dest = 'snps_vcf', required = True)
    _snp_pileup.add_argument('--normal-bam', dest = 'normal_bam', required = True)
    _snp_pileup.add_argument('--tumor-bam', dest = 'tumor_bam', required = True)
    _snp_pileup.add_argument('--output-prefix', dest = 'output_prefix', required = True)
    _snp_pileup.set_defaults(func = SnpPileup._run)
    """
    $ ./run.py snp_pileup \
    --tumor-bam /juno/work/ci/helix_filters_01/fixtures/Proj_08390_G/bam/Sample34.rg.md.abra.printreads.bam \
    --normal-bam /juno/work/ci/helix_filters_01/fixtures/Proj_08390_G/bam/Sample33.rg.md.abra.printreads.bam \
    --snps-vcf /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf \
    --output-prefix Sample34.Sample33
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
