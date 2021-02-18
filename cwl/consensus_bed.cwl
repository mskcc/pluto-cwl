#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: "
Workflow to merge a large number of maf files into a single consensus bed file
"
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  maf_files:
    type: File[]

steps:
  convert_mafs:
    run: maf2bed.cwl
    scatter: file
    in:
      file: maf_files
      maf_file:
        valueFrom: ${ return inputs.file; }
    out:
      [ output_file ]

  merge_beds:
    run: mergebed.cwl
    in:
      bed_files: convert_mafs/output_file
    out:
      [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: merge_beds/output_file
