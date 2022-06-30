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
  msi_files:
    type: File[]
    doc: "msi.tsv files"


steps:  

  # concatenate all the individual MSI tables into a single table
  concat_msi_tables:
    run: concat-tables.cwl
    in:
      input_files: msi_files
      output_filename:
        valueFrom: ${ return "msi.tsv"; }
      comments:
        valueFrom: ${ return true; }
    out:
      [ output_file ]

  #replace % colname with MSI_SCORE
  replace_col_name:
    run: replace_colname.cwl
    in:
      old_name:
        valueFrom: ${ return "%"; }
      new_name:
        valueFrom: ${ return "MSI_SCORE"; }
      input_file: concat_msi_tables/output_file
      output_filename:
        valueFrom: ${ return "msi-replaced.tsv"; }
    out:
      [ output_file ]

  # add msi status based on the msi score
  add_msi_status_col:
    run: add_msi_status.cwl
    in:
      input_filename: replace_col_name/output_file
      output_filename:
        valueFrom: ${ return "msi_status_added.tsv"; }
      header:
        valueFrom: ${ return "MSI_STATUS"; }
    out:
      [ output_file ]


  # cut msi.tsv to get only msi scores
  cut_msi_table:
    run: cut.cwl
    in:
      field_indexes:
        valueFrom: ${ return "3,4,5"; }
      input_file: add_msi_status_col/output_file
    out:
      [ output_file ]


  # combine the MSI results with the data clinical file
  merge_data_clinical:
    run: merge-tables.cwl
    in:
      table1: data_clinical_file
      table2: cut_msi_table/output_file
      key1:
        valueFrom: ${ return "SAMPLE_ID"; } # sample column header from data clinical file
      key2:
        valueFrom: ${ return "SAMPLE_ID"; } # sample column header from MSI file
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
