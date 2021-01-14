#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: "
Workflow for calculating TMB tumor mutational burden
"
requirements:
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  mutations_file:
    type: File
    doc: "File with mutations for the sample"
  assay_coverage:
    type: string
    doc: "genome_coverage value; amount of the genome in bp covered by the assay"
  sample_id:
    type: string

steps:
    filter_variants:
      run: tmb_variant_filter.cwl
      in:
        input_file: mutations_file
        output_filename:
          valueFrom: ${ return "filter_variants.maf"; }
      out:
        [ output_file ]

    calc_tmb_value:
      run: calc-tmb.cwl
      in:
        input_file: filter_variants/output_file
        output_filename:
          valueFrom: ${ return "tmb.txt"; }
        genome_coverage: assay_coverage
      out:
        [ output_file ]

    fix_tmb_header:
      run: add_header.cwl
      in:
        input_file: calc_tmb_value/output_file
        header_str:
          valueFrom: ${ return "TMB"; }
      out:
        [ output_file ]

    add_sampleID:
      run: paste-col.cwl
      in:
        input_file: fix_tmb_header/output_file
        output_filename:
          valueFrom: ${ return "tmb.tsv"; }
        header:
          valueFrom: ${ return "SampleID"; }
        value: sample_id
      out:
        [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: add_sampleID/output_file
