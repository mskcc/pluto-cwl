#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
doc: "
Example CWL workflow that uses some advanced features
"
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  value:
    type: string
    doc: "a value to pass in everywhere"
  samples:
    type:
      type: array
      items:
        type: record
        fields:
          sample_id: string

steps:
  # run example.cwl for each sample
  run_example:
    run: example.cwl
    scatter: sample
    in:
      sample: samples
      sample_id:
        valueFrom: ${ return inputs.sample['sample_id']; }
      value: value
    out:
      [ output_file ] # [ output.tsv 1, output.tsv 2, ... ] array of output table files for each input sample

  # concatenate all the individual output tables into a single table
  concat_tables:
    run: concat-tables.cwl
    in:
      input_files: run_example/output_file # array of input files
      output_filename:
        valueFrom: ${ return "output.concat.tsv"; }
    out:
      [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: concat_tables/output_file
