#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: "
Wrapper to run indexing on all bams before submitting for samples fillout
Includes secondary input channels to allow for including .bam files that do not have indexes

NOTE: need v1.1 upgrade so we can do it all from a single channel with optional secondary files;
https://www.commonwl.org/v1.1/CommandLineTool.html#SecondaryFileSchema
"
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  sample_ids:
    type:
        type: array
        items: string

  bam_files:
    type:
        type: array
        items: File
    secondaryFiles:
        - ^.bai

  maf_files:
    type:
      type: array
      items: File

  unindexed_bam_files:
    type:
      type: array
      items: File

  unindexed_sample_ids:
    type:
      type: array
      items: string

  unindexed_maf_files:
    type:
      type: array
      items: File

  ref_fasta:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
      - .fai
      - ^.dict

  exac_filter: # need this to resolve error in subworkflow: Anonymous file object must have 'contents' and 'basename' fields.
    type: File
    default:
      class: File
      path: /juno/work/ci/resources/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz

  is_impact:
    type: boolean
    default: True

  argos_version_string:
    type: [ "null", string ]
    default: "Unspecified"

  fillout_output_fname:
    type: string
    default: "fillout.maf"

steps:
  # index files
  run_indexer:
    run: index_bam.cwl
    in:
      bam: unindexed_bam_files
    scatter: [ bam ]
    scatterMethod: dotproduct
    out: [ bam_indexed ]

  run_maf_filter:
    run: maf_filter.cwl
    in:
      maf_file: maf_files
      is_impact: is_impact
      argos_version_string: argos_version_string
    scatter: [ maf_file ]
    scatterMethod: dotproduct
    out: [ analysis_mutations_file ]

  run_samples_fillout:
    run: samples_fillout_workflow.cwl
    in:
      output_fname: fillout_output_fname
      exac_filter: exac_filter
      sample_ids:
        source: [ sample_ids, unindexed_sample_ids ]
        linkMerge: merge_flattened
      bam_files:
        source: [ bam_files, run_indexer/bam_indexed ]
        linkMerge: merge_flattened
      maf_files:
        source: [ run_maf_filter/analysis_mutations_file, unindexed_maf_files ]
        linkMerge: merge_flattened
      ref_fasta: ref_fasta
    out: [ output_file ]

outputs:

  output_file:
    type: File
    outputSource: run_samples_fillout/output_file
