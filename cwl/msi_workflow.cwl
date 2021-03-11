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
          normal_bam:
            type: File
            secondaryFiles:
              - ^.bai
          tumor_bam:
            type: File
            secondaryFiles:
              - ^.bai


steps:
  run_msi_add_sample_id:
    in:
      microsatellites_file: microsatellites_file
      pair: pairs
      normal_bam:
        valueFrom: ${ return inputs.pair['normal_bam']; }

      tumor_bam:
        valueFrom: ${ return inputs.pair['tumor_bam']; }

      col_header:
        valueFrom: ${ return 'SAMPLE_ID'; }
      tumor_id:
        valueFrom: ${ return inputs.pair['tumor_id']; }

    scatter: pair

    out: [ output_file ]


    run:
      class: Workflow
      id: run_msi_and_add_sample_id
      inputs:
        microsatellites_file: File
        normal_bam: File
        tumor_bam: File
        tumor_id: string

      outputs:
        output_file:
          type: File
          outputSource: add_sample_id/output_file

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

  # concatenate all the individual MSI tables into a single table
  concat_msi_tables:
    run: concat-tables.cwl
    in:
      input_files: run_msi_add_sample_id/output_file
      output_filename:
        valueFrom: ${ return "msi.tsv"; }
      comments:
        valueFrom: ${ return true; }
    out:
      [ output_file ]

  # cut msi.tsv to get only msi scores
  cut_msi_table:
    run: cut.cwl
    in:
      field_indexes:
        valueFrom: ${ return "3,4"; }
      input_file: concat_msi_tables/output_file
    out:
      [ output_file ]

  #replace % colname with MSI_SCORE
  replace_col_name:
    run: replace2.cwl
    in:
      input_file: cut_msi_table/output_file
      output_filename:
        valueFrom: ${ return "msi-replaced.tsv"; }
    out:
      [ output_file ]


  # combine the MSI results with the data clinical file
  merge_data_clinical:
    run: merge-tables.cwl
    in:
      table1: data_clinical_file
      table2: replace_col_name/output_file
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
