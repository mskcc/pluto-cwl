#!/usr/bin/env cwl-runner

cwlVersion: v1.0
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

# input for
# cwlVersion: v1.1
# NOTE: Important! Need  cwlVersion: v1.1 for the array record fields secondaryFiles to work here
# inputs:
#   maf_file:
#     type: File
#     doc: "The reference maf file with coordinates to run fillout against for each bam file"
#   bams:
#     type:
#       type: array
#       items:
#         type: record
#         fields:
#           bam_file:
#             type: File
#             secondaryFiles:
#               - ^.bai
#           sample_id: string
#   ref_fasta:
#     type: File
#     secondaryFiles:
#       - .amb
#       - .ann
#       - .bwt
#       - .pac
#       - .sa
#       - .fai
#       - ^.dict

inputs:
  maf_file:
    type: File
    doc: "The reference maf file with coordinates to run fillout against for each bam file"
  sample_ids:
    type:
        type: array
        items: string
  bam_files:
    type:
        type: array
        items: File
    secondaryFiles:
        - ^.bai
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
    scatter: [ sample_id, bam ]
    scatterMethod: dotproduct
    in:
      bam: bam_files
      sample_id: sample_ids
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
