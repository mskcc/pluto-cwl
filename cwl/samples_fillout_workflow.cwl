#!/usr/bin/env cwl-runner

# NOTE: A newer version that will deprecate this exists once cwlVersion gets updated to 1.1
cwlVersion: v1.0
class: Workflow
doc: "
Workflow to run GetBaseCountsMultiSample fillout on a number of samples, each with their own bam and maf files
"
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  # NOTE: arrays for sample_ids, bam_files, maf_files must all be the same length and in the same order by sample
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
  maf_files:
    type:
        type: array
        items: File
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


# # future version of input format;
# NOTE: Important! Need  cwlVersion: v1.1 for the array record fields secondaryFiles to work here
# cwlVersion: v1.1
# inputs:
#   samples:
#     type:
#       type: array
#       items:
#         type: record
#         fields:
#           bam_file:
#             type: File
#             secondaryFiles:
#               - ^.bai
#           maf_file: File
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



steps:

  # combine all the sample maf files into a single reference maf to run against GetBaseCountsMultiSample
  make_reference_maf:
    run: concat-mafs.cwl
    in:
      input_files: maf_files
    out:
      [ output_file ]

  # check the base counts at each position in the sample reference maf for each sample
  sample_fillout:
    run: getbasecountsmultisample.cwl
    scatter: [ sample_id, bam ]
    scatterMethod: dotproduct
    in:
      sample_id: sample_ids
      bam: bam_files
      ref_fasta: ref_fasta
      maf: make_reference_maf/output_file
    out:
      [ output_file ] # [ "fillout.maf", ... ]

  # merge all the samples maf fillouts
  combine_fillouts:
    run: concat-mafs_all_cols.cwl
    in:
      input_files: sample_fillout/output_file
    out:
      [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: combine_fillouts/output_file
