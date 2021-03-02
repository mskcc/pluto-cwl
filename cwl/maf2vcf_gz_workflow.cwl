#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: "
Workflow to convert a maf file into a vcf.gz with .tbi index
"
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  maf_file: File
  ref_fasta:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
      - .fai
      - ^.dict

steps:
  maf2vcf:
    run: maf2vcf.cwl
    in:
      maf_file: maf_file
      ref_fasta: ref_fasta
    out:
      [ output_file ]

  sort_vcf:
    run: vcf_sort.cwl
    in:
      vcf_file: maf2vcf/output_file
    out:
      [ output_file ]

  zip_vcf:
    run: bgzip.cwl
    in:
      input_file: sort_vcf/output_file
      output_filename:
        valueFrom: ${ return "variants.vcf.gz"; }
    out:
      [ output_file ]

  index_vcf:
    run: index_vcf.cwl
    in:
      input_file: zip_vcf/output_file
    out:
      [ indexed_file ]

outputs:
  output_file:
    type: File
    outputSource: index_vcf/indexed_file
