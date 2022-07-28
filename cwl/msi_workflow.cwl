#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
doc: "
Workflow to run the MSI analysis on a batch of samples
"
requirements:
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
  - $import: types.yml

inputs:
  threads:
    type: string
    default: "8"
  microsatellites_file:
    doc: reference microsatellite list
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
  # NOTE: these two arrays of File with secondaryFiles should eventually be merged directly into the `pairs` record array
  # after upgrading Toil to support cwlVersion 1.1
  normal_bam_files:
    type:
        type: array
        items: File
    secondaryFiles:
        - ^.bai

  tumor_bam_files:
    type:
        type: array
        items: File
    secondaryFiles:
        - ^.bai


steps:
  run_msi_add_sample_id:
    scatter: [ pair, normal_bam, tumor_bam ]
    scatterMethod: dotproduct
    in:
      threads: threads
      microsatellites_file: microsatellites_file
      pair: pairs
      normal_bam: normal_bam_files
      tumor_bam: tumor_bam_files
      pair_id:
        valueFrom: ${ return inputs.pair['pair_id']; }
      tumor_id:
        valueFrom: ${ return inputs.pair['tumor_id']; }
      normal_id:
        valueFrom: ${ return inputs.pair['normal_id']; }
    out: [ pair ]
    run:
      class: Workflow
      inputs:
        threads: string
        microsatellites_file: File
        normal_bam: File
        tumor_bam: File
        pair_id: string
        tumor_id: string
        normal_id: string
      outputs:
        pair:
          type: "types.yml#MSIOutputPair"
          outputSource: create_msi_pair_output/pair
        # output_file:
        #   type: File
        #   outputSource: add_sample_id/output_file
      steps:
        # run the MSI analysis for each tumor sample in the list of pairs
        run_msi_workflow:
          doc: run msisensor for microsatellite instability analysis
          run: msi.cwl
          in:
            threads: threads
            microsatellites_file: microsatellites_file
            normal_bam: normal_bam
            tumor_bam: tumor_bam
          out:
            [output_file]
        # output_file looks like this;
        # Total_Number_of_Sites   Number_of_Somatic_Sites %
        # 628     138     21.97

        replace_col_name:
          doc: replace colname with MSI_SCORE
          run: replace_colname.cwl
          in:
            old_name:
              valueFrom: ${ return "%"; }
            new_name:
              valueFrom: ${ return "MSI_SCORE"; }
            input_file: run_msi_workflow/output_file
            output_filename:
              valueFrom: ${ return "msi-replaced.tsv"; }
          out:
            [ output_file ]

        add_msi_status_col:
          doc: add msi status based on the msi score
          run: add_msi_status.cwl
          in:
            input_filename: replace_col_name/output_file
            output_filename:
              valueFrom: ${ return "msi_status_added.tsv"; }
            header:
              valueFrom: ${ return "MSI_STATUS"; }
          out:
            [ output_file ]

        cut_msi_table:
          doc: cut msi.tsv to get only msi scores
          run: cut.cwl
          in:
            field_indexes:
              valueFrom: ${ return "3,4,5"; }
            input_file: add_msi_status_col/output_file
          out:
            [ output_file ]

        add_sample_id:
          doc: add a column for sample ID to the file
          run: paste-col.cwl
          in:
            pair_id: pair_id
            input_file: cut_msi_table/output_file
            output_filename:
              valueFrom: ${ return inputs.pair_id + '.msi.tsv'; }
            header:
              valueFrom: ${ return 'SAMPLE_ID'; }
            value: tumor_id
          out:
            [ output_file ]
          # output_file looks like this;
          # MSI_SCORE       MSI_STATUS      SAMPLE_ID
          # 21.97   Instable        Sample1-T

        create_msi_pair_output:
          doc:
          in:
            msi_tsv: add_sample_id/output_file
            pair_id: pair_id
            tumor_id: tumor_id
            normal_id: normal_id
          out: [ pair ]
          run:
            class: ExpressionTool
            inputs:
              msi_tsv: File
              pair_id: string
              tumor_id: string
              normal_id: string
            outputs:
              pair: "types.yml#MSIOutputPair"
            expression: |
              ${
                var pair = {
                  "msi_tsv": inputs.msi_tsv,
                  "pair_id": inputs.pair_id,
                  "tumor_id": inputs.tumor_id,
                  "normal_id": inputs.normal_id,
                };
                return {
                  "pair": pair
                };
              }


outputs:
  pairs:
    type: "types.yml#MSIOutputPair[]"
    outputSource: run_msi_add_sample_id/pair
