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
    run: msi.cwl
    scatter: pair
    in:
      pair: pairs
      mutations_file:
        valueFrom: ${ return inputs.pair['pair_maf']; }
      sample_id:
        valueFrom: ${ return inputs.pair['tumor_id']; }
      assay_coverage: assay_coverage
    out:
      [ output_file ]

  # # concatenate all the individual TMB tables into a single table
  # concat_tmb_tables:
  #   run: concat-tables.cwl
  #   in:
  #     input_files: run_tmb_workflow/output_file
  #     output_filename:
  #       valueFrom: ${ return "tmb.tsv"; }
  #     comments:
  #       valueFrom: ${ return true; }
  #   out:
  #     [ output_file ]
  #
  # # combine the TMB results with the data clinical file
  # merge_data_clinical:
  #   run: merge-tables.cwl
  #   in:
  #     table1: data_clinical_file
  #     table2: concat_tmb_tables/output_file
  #     key1:
  #       valueFrom: ${ return "SAMPLE_ID"; } # sample column header from data clinical file
  #     key2:
  #       valueFrom: ${ return "SampleID"; } # sample column header from TMB file
  #     output_filename:
  #       valueFrom: ${ return "data_clinical_sample.txt"; } # TODO: should this be passed in?
  #   out:
  #     [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: merge_data_clinical/output_file
