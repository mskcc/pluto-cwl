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
  microsatellites_file: File
  normal_bam: File
  tumor_bam: File
  tumor_id: string

steps:
  # run the MSI analysis for each tumor sample in the list of pairs
  run_msi_workflow:
    run: msi.cwl
    in:
      d: microsatellites_file
      n: normal_bam
      t: tumor_bam
      o:
        valueFrom: ${ return 'msi.txt'; }
    out:
      [output_file]

  add_sample_id:
    run: paste-col.cwl
    in:
      input_file: run_msi_workflow/output_file
      output_filename:
        valueFrom: ${ return 'msi.tsv'; }
      header:
        valueFrom: ${ return 'SAMPLE_ID'; }
      value: tumor_id
    out:
      [output_file]

outputs:
  output_file:
    type: File
    outputSource: add_sample_id/output_file
