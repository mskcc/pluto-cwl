#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
doc: "
Workflow to run the TMB analysis on a batch of samples and merge the results back into a single data clinical file
"
requirements:
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
  - $import: types.yml

inputs:
  assay_coverage:
    type: string
    doc: "genome_coverage value; amount of the genome in bp covered by the assay"
  pairs:
    type: "types.yml#TMBInputPair[]"

steps:
  # run the TMB analysis for each tumor sample in the list of pairs
  run_tmb_workflow:
    scatter: pair
    in:
      pair: pairs
      mutations_file:
        valueFrom: ${ return inputs.pair['pair_maf']; }
      sample_id:
        valueFrom: ${ return inputs.pair['tumor_id']; }
      normal_id:
        valueFrom: ${ return inputs.pair['normal_id']; }
      pair_id:
        valueFrom: ${ return inputs.pair['pair_id']; }
      assay_coverage: assay_coverage
    out:
      [ pair ]
    run:
      class: Workflow
      inputs:
        mutations_file:
          type: File
          doc: "File with mutations for the sample"
        assay_coverage:
          type: string
          doc: "genome_coverage value; amount of the genome in bp covered by the assay"
        sample_id:
          type: string
        normal_id:
          type: string
        pair_id:
          type: string
      outputs:
        pair:
          type: "types.yml#TMBOutputPair"
          outputSource: create_tmb_pair_output/pair
      steps:
        filter_variants:
          doc: filter the variant maf file for only the variants desired for use in TMB calculation
          run: tmb_variant_filter.cwl
          in:
            pair_id: pair_id
            input_file: mutations_file
            output_filename:
              valueFrom: ${ return inputs.pair_id + ".tmb.maf"; }
          out:
            [ output_file ]

        calc_tmb_value:
          doc: calculate the TMB for the variants present based on assay coverage
          run: calc-tmb.cwl
          in:
            pair_id: pair_id
            input_file: filter_variants/output_file
            output_filename:
              valueFrom: ${ return inputs.pair_id + ".tmb.txt"; }
            genome_coverage: assay_coverage
            normal_id: normal_id
          out:
            [ output_file ]

        fix_tmb_header:
          doc: turn the TMB value into a table format with header
          run: add_header.cwl
          in:
            input_file: calc_tmb_value/output_file
            header_str:
              valueFrom: ${ return "CMO_TMB_SCORE"; }
          out:
            [ output_file ]

        add_sampleID:
          doc: add the sample ID back to the TMB table file
          run: paste-col.cwl
          in:
            pair_id: pair_id
            input_file: fix_tmb_header/output_file
            output_filename: # NOTE: we plan to concat this file later so it needs to have a unique filename !!
              valueFrom: ${ return inputs.pair_id + ".tmb.tsv"; }
            header:
              valueFrom: ${ return "SampleID"; }
            value: sample_id
          out:
            [ output_file ]

        add_assay_coverage:
          doc: add the assay coverage to the table
          run: paste-col.cwl
          in:
            pair_id: pair_id
            input_file: add_sampleID/output_file
            output_filename: # NOTE: we plan to concat this file later so it needs to have a unique filename !!
              valueFrom: ${ return inputs.pair_id + ".tmb.tsv"; }
            header:
              valueFrom: ${ return "CMO_ASSAY_COVERAGE"; }
            value: assay_coverage
          out:
            [ output_file ]

        create_tmb_pair_output:
          doc: gather the TMB analysis outputs into a pair entry
          in:
            pair_id: pair_id
            tumor_id: sample_id
            normal_id: normal_id
            tmb_maf: filter_variants/output_file
            tmb_tsv: add_assay_coverage/output_file
          out: [ pair ]
          run:
            class: ExpressionTool
            inputs:
              pair_id: string
              tumor_id: string
              normal_id: string
              tmb_maf: File
              tmb_tsv: File
            outputs:
              pair: "types.yml#TMBOutputPair"
            expression: |
              ${
                var pair = {
                  "pair_id": inputs.pair_id,
                  "tumor_id": inputs.tumor_id,
                  "normal_id": inputs.normal_id,
                  "tmb_maf": inputs.tmb_maf,
                  "tmb_tsv": inputs.tmb_tsv,
                };

                return {"pair": pair};
              }


outputs:
  pairs:
    type: "types.yml#TMBOutputPair[]"
    outputSource: run_tmb_workflow/pair
