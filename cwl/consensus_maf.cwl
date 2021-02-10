#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: "
Workflow to merge a large number of maf files into a single consensus maf file for use with GetBaseCountsMultiSample
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
  # concatenate all the mafs into a single file for use with GetBaseCountsMultiSample
  merge_all_mafs:
    run: concat-mafs.cwl
    in:
      input_files: maf_files
    out:
      [ output_file ] # output.maf

  # deduplicate the maf entries, and sort
  deduplicate_maf:
    run: deduplicate-maf.cwl
    in:
      input_file: merge_all_mafs/output_file
    out:
      [ output_file ] # dedup.maf

outputs:
  output_file:
    type: File
    outputSource: deduplicate_maf/output_file
