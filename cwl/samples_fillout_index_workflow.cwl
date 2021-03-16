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

steps:

  # index files
  run_indexer:
    run: index_bam.cwl
    in:
      bam: unindexed_bam_files
    scatter: [ bam ]
    scatterMethod: dotproduct
    out: [ bam_indexed ]

  run_samples_fillout:
    run: samples_fillout_workflow.cwl
    in:
      sample_ids:
        source: [ sample_ids, unindexed_sample_ids ]
        linkMerge: merge_flattened
      bam_files:
        source: [ bam_files, run_indexer/bam_indexed ]
        linkMerge: merge_flattened
      maf_files:
        source: [ maf_files, unindexed_maf_files ]
        linkMerge: merge_flattened
      ref_fasta: ref_fasta
    out: [ output_file ]

outputs:

  output_file:
    type: File
    outputSource: run_samples_fillout/output_file
