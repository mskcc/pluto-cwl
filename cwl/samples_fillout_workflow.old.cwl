#!/usr/bin/env cwl-runner

# NOTE: Important! Need  cwlVersion: v1.1 for the array record fields secondaryFiles to work here
cwlVersion: v1.1
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
  samples:
    type:
      type: array
      items:
        type: record
        fields:
          bam_file:
            type: File
            secondaryFiles:
              - ^.bai
          maf_file: File
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
  # need to get all the maf files from all the input samples
  get_all_mafs:
    in:
      samples: samples
    out:
      [ maf_files ]
    run:
      class: ExpressionTool
      id: get_all_mafs
      inputs:
        samples:
          type:
            type: array
            items:
              type: record
              fields:
                maf_file: File
      outputs:
        maf_files: File[]
      expression: "${
        var maf_file_list = [];
        for(var i in inputs.samples){
          maf_file_list.push(inputs.samples[i]['maf_file'])
        };
        return {'maf_files':maf_file_list};
      }"

  # combine all the sample maf files into a single reference maf to run against GetBaseCountsMultiSample
  make_reference_maf:
    run: concat-mafs.cwl
    in:
      input_files: get_all_mafs/maf_files
    out:
      [ output_file ]

  # check the base counts at each position in the sample reference maf for each sample
  sample_fillout:
    run: getbasecountsmultisample.cwl
    scatter: sample
    in:
      sample: samples
      bam:
        valueFrom: ${ return inputs.sample['bam_file']; }
      sample_id:
        valueFrom: ${ return inputs.sample['sample_id']; }
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
