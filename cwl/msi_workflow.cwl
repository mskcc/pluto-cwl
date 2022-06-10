#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: "
Workflow to run the MSI analysis on a batch of samples and merge the results back into a single data clinical file
"
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  data_clinical_file:
    type: File
    doc: "data clinical samplesheet file to merge the MSI results into"
  microsatellites_file:
    type: File
  pairs:
    type:
      type: array
      items:
        type: record
        fields:
          pair_id: string
          tumor_id: string
          normal_id: string
  # NOTE: these two arrays of File with secondaryFiles should eventually be merged directly into the `pairs` record array
  # after upgrading Toil to support cwlVersion 1.1
  normal_bam_files:
    type:
        type: array
        items: File
    secondaryFiles:
        - ^.bai

  tumor_bam_files:
    type:
        type: array
        items: File
    secondaryFiles:
        - ^.bai


steps:
  run_msi_add_sample_id:
    scatter: [ pair, normal_bam, tumor_bam ]
    scatterMethod: dotproduct
    run: run_msi_and_add_sample_id.cwl
    in:
      microsatellites_file: microsatellites_file
      pair: pairs
      normal_bam: normal_bam_files
      tumor_bam: tumor_bam_files
      tumor_id:
        valueFrom: ${ return inputs.pair['tumor_id']; }
    out: [ output_file ]
  
  merge_msi:
    run: merge_msi.cwl
    in:
      data_clinical_file: data_clinical_file
      msi_files: run_msi_add_sample_id/output_file
    out: [ output_file ] 

outputs:
  output_file:
    type: File
    outputSource: merge_msi/output_file
