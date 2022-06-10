#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: "
Workflow to run the TMB analysis on a batch of samples and merge the results back into a single data clinical file
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
    doc: "data clinical samplesheet file to merge the TMB results into"
  assay_coverage:
    type: string
    doc: "genome_coverage value; amount of the genome in bp covered by the assay"
  pairs:
    type:
      type: array
      items:
        type: record
        fields:
          pair_maf: File
          pair_id: string
          tumor_id: string
          normal_id: string

steps:
  # run the TMB analysis for each tumor sample in the list of pairs
  run_tmb_workflow:
    run: tmb.cwl
    scatter: pair
    in:
      pair: pairs
      mutations_file:
        valueFrom: ${ return inputs.pair['pair_maf']; }
      sample_id:
        valueFrom: ${ return inputs.pair['tumor_id']; }
      normal_id:
        valueFrom: ${ return inputs.pair['normal_id']; }
      assay_coverage: assay_coverage
    out:
      [ output_file ] # [ tmb.tsv 1, tmb.tsv 2, ... ] array of tmb table files for each input pair
  
  merge_tmb:
    run: merge_tmb.cwl
    in:
      data_clinical_file: data_clinical_file
      tmb_files: run_tmb_workflow/output_file
    out:
      [ output_file ]
outputs:
  output_file:
    type: File
    outputSource: merge_tmb/output_file
