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
  tmb_files:
    type: File[]
    doc: "*.tmb.tsv files"

steps:

  # concatenate all the individual TMB tables into a single table
  concat_tmb_tables:
    run: concat-tables_dir.cwl # NOTE: Important!! use this CWL in case a huge amount of files are passed!!
    in:
      input_files: tmb_files
      output_filename:
        valueFrom: ${ return "tmb.concat.tsv"; }
      comments:
        valueFrom: ${ return true; }
    out:
      [ output_file ]

  # combine the TMB results with the data clinical file
  merge_data_clinical:
    run: merge-tables.cwl
    in:
      table1: data_clinical_file
      table2: concat_tmb_tables/output_file
      key1:
        valueFrom: ${ return "SAMPLE_ID"; } # sample column header from data clinical file
      key2:
        valueFrom: ${ return "SampleID"; } # sample column header from TMB file
      output_filename:
        valueFrom: ${ return "data_clinical_sample.txt"; } # TODO: should this be passed in?
      cBioPortal:
        valueFrom: ${ return true; }
    out:
      [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: merge_data_clinical/output_file
