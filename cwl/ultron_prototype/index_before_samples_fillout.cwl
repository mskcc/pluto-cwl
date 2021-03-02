#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: "
Wrapper to run indexing on all bams before submitting for samples fillout
"
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  sample_names:
    type:
        type: array
        items: string

  bam_files:
    type:
        type: array
        items: File
    secondaryFiles:
        - ^.bai
        - .bai

  dmp_bams:
    type:
      type: array
      items: File

  dmp_bams_sample_names:
    type:
      type: array
      items: string

  maf_files: File[]

  dmp_maf_files: File[]

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
    run: cmo_index.cwl 
    in:
      bam: dmp_bams
    scatter: [ bam ]
    scatterMethod: dotproduct
    out: [ bam_indexed ]

  run_samples_fillout:
    run: samples_fillout_workflow_ultron.cwl
    in:
      sample_names: sample_names
      dmp_bams_sample_names: dmp_bams_sample_names
      bam_files: bam_files
      maf_files: maf_files
      dmp_maf_files: dmp_maf_files
      sample_name:
        source: [ sample_names, dmp_bams_sample_names ]
        linkMerge: merge_flattened
      bam_file:
        source: [ bam_files, run_indexer/bam_indexed ]
        linkMerge: merge_flattened
      maf_file: 
        source: [ maf_files, dmp_maf_files ]
        linkMerge: merge_flattened
      ref_fasta: ref_fasta
    out: [ output_file ]

outputs:

  output_file:
    type: File
    outputSource: run_samples_fillout/output_file
