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
          pair_normal_bam:
            type: File
            secondaryFiles:
              - ^.bai
          pair_tumor_bam:
            type: File
            secondaryFiles:
              - ^.bai

  normal_sample_name:
    valueFrom: ${ return inputs.pair['pair_normal_bam']; }
  tumor_sample_name:
    valueFrom: ${ return inputs.pair['pair_tumor_bam']; }
  output_file_name:
    valueFrom: ${ return inputs.pair['pair_id']; }
  col_header:
    valueFrom: ${ return 'SAMPLE_ID'; }
  col_value:
    valueFrom: ${ return inputs.pair['pair_id']; }

outputs:
  output_file:
    type: File
    outputSource: concat_msi_tables/output_file # merge_data_clinical/output_file


steps:
  run_msi_add_sample_id:
    in:
      pair: pairs
      normal_sample_name: normal_sample_name
      tumor_sample_name: tumor_sample_name
      output_file_name: output_file_name
      col_header: col_header
      col_value: col_value
      # pair: pairs
      # normal_sample_name: File
      # tumor_sample_name: File
      # output_file_name: string
      # col_header: string
      # col_value: string

    scatter: pair

    out: [ output_file ]


    run:
      class: Workflow
      id: run_msi_and_add_sample_id
      inputs:
        normal_sample_name: File
        tumor_sample_name: File
        output_file_name: string
        col_header: string
        col_value: string
        # pair: pairs
        # normal_sample_name: valueFrom: ${ return inputs.pair['pair_normal_bam']; }
        # tumor_sample_name: valueFrom: ${ return inputs.pair['pair_tumor_bam']; }
        # output_file_name: valueFrom: ${ return inputs.pair['pair_id']; }
        # col_header: valueFrom: ${ return 'SAMPLE_ID'; }
        # col_value: valueFrom: ${ return inputs.pair['pair_id']; }

      outputs: add_sample_id/output_file
      steps:
        # run the MSI analysis for each tumor sample in the list of pairs
        run_msi_workflow:
          run: msi.cwl
          in:
            d: microsatellites_file
            n: normal_sample_name
            t: tumor_sample_name
            o: output_file_name
          out:
            [output_file]

        add_sample_id:
          run: paste-col.cwl
          in:
            input_file: run_msi_workflow/output_file
            output_filename: output_file_name
            header: col_header
            value: col_value
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

  # combine the MSI results with the data clinical file
  # merge_data_clinical:
  #   run: merge-tables.cwl
  #   in:
  #     table1: data_clinical_file
  #     table2: concat_msi_tables/output_file
  #     key1:
  #       valueFrom: ${ return "SAMPLE_ID"; } # sample column header from data clinical file
  #     key2:
  #       valueFrom: ${ return "SampleID"; } # sample column header from MSI file
  #     output_filename:
  #       valueFrom: ${ return "data_clinical_sample.txt"; } # TODO: should this be passed in?
  #   out:
  #     [ output_file ]
