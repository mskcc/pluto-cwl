#!/usr/bin/env cwl-runner

cwlVersion: v1.1
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
          pair_id: string
          tumor_id: string
          normal_id: string
          pair_normal_bam:
            type: File
            secondaryFiles:
              - ^.bai
          pair_tumor_bam:
            type: File
            secondaryFiles:
              - ^.bai

steps:
  # run the MSI analysis for each tumor sample in the list of pairs
  run_msi_workflow:
    run: msi.cwl
    scatter: pair
    in:
      pair: pairs
      n:
        valueFrom: ${ return inputs.pair['pair_normal_bam']; }
      t:
        valueFrom: ${ return inputs.pair['pair_tumor_bam']; }
      d:
        valueFrom: ${ return '/work/ci/resources/request_files/msisensor/b37_known_somatic_microsatellites.list'; }
      o:
        valueFrom: ${ return inputs.pair['pair_id']; }
    out:
      [ output_file ]

  # concatenate all the individual MSI tables into a single table
  concat_msi_tables:
    run: concat-tables.cwl
    in:
      input_files: run_msi_workflow/output_file
      output_filename:
        valueFrom: ${ return "msi.tsv"; }
      comments:
        valueFrom: ${ return true; }
    out:
      [ output_file ]

  # combine the MSI results with the data clinical file
  merge_data_clinical:
    run: merge-tables.cwl
    in:
      table1: data_clinical_file
      table2: concat_msi_tables/output_file
      key1:
        valueFrom: ${ return "SAMPLE_ID"; } # sample column header from data clinical file
      key2:
        valueFrom: ${ return "SampleID"; } # sample column header from MSI file
      output_filename:
        valueFrom: ${ return "data_clinical_sample.txt"; } # TODO: should this be passed in?
    out:
      [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: merge_data_clinical/output_file
