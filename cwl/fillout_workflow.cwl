#!/usr/bin/env cwl-runner

# NOTE: Important! Need  cwlVersion: v1.1 for the array record fields secondaryFiles to work here
cwlVersion: v1.1
class: Workflow
doc: "
Workflow to run GetBaseCountsMultiSample fillout on a number of bam files with a single maf file
"
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  maf_file:
    type: File
    doc: "The reference maf file with coordinates to run fillout against for each bam file"
  bams:
    type:
      type: array
      items:
        type: record
        fields:
          bam_file:
            type: File
            secondaryFiles:
              - ^.bai
          sample_id: string
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
  fillout:
    run: getbasecountsmultisample.cwl
    scatter: bam_record
    in:
      bam_record: bams
      bam:
        valueFrom: ${ return inputs.bam_record['bam_file']; }
      sample_id:
        valueFrom: ${ return inputs.bam_record['sample_id']; }
      ref_fasta: ref_fasta
      maf: maf_file
    out:
      [ output_file ] # [ "fillout.maf", ... ]

  concat_mafs:
    run: concat-mafs_all_cols.cwl
    in:
      input_files: fillout/output_file
    out:
      [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: concat_mafs/output_file
