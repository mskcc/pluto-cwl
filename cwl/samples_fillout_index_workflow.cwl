#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: "
Wrapper to run indexing on all bams before submitting for samples fillout
Includes secondary input channels to allow for including .bam files that do not have indexes
Also include other extra handling needed for files that might not meet needs for the fillout workflow

NOTE: need v1.1 upgrade so we can do it all from a single channel with optional secondary files;
https://www.commonwl.org/v1.1/CommandLineTool.html#SecondaryFileSchema
"
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  samples: # NOTE: in prod, these end up being the research samples
    type:
      type: array
      items:
        type: record
        fields:
          maf_file: File
          sample_id: string # must match sample ID used inside maf file
          normal_id: string
  bam_files:
    type:
        type: array
        items: File
    secondaryFiles:
        - ^.bai

  unindexed_samples: # NOTE: in prod, these end up being the clinical samples
    type:
      type: array
      items:
        type: record
        fields:
          maf_file: File
          sample_id: string # must match sample ID used inside maf file
          normal_id: string

  unindexed_bam_files:
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
  exac_filter: # need this to resolve error in subworkflow: Anonymous file object must have 'contents' and 'basename' fields.
  # TODO: this needs the .tbi/.csi index file added!!
    type: File
    default:
      class: File
      path: /juno/work/ci/resources/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz

  # these are needed for the filter script
  is_impact:
    type: boolean
    default: True

  argos_version_string:
    type: [ "null", string ]
    default: "Unspecified"

  fillout_output_fname:
    type: string
    default: "fillout.maf"

steps:
  # index any bam files that lacked .bai indexes
  run_indexer:
    run: index_bam.cwl
    in:
      bam: unindexed_bam_files
    scatter: bam
    out: [ bam_indexed ]

  # run filter script to apply cBioPortal filters on variants for fillout
  run_maf_filter:
    run: maf_filter.cwl
    in:
      sample: samples
      maf_file:
        valueFrom: ${ return inputs.sample['maf_file']; }
      is_impact: is_impact
      argos_version_string: argos_version_string
    scatter: sample
    out: [ cbio_mutation_data_file ]

  # NOTE: In prod, the unindexed_samples end up being the clinical samples; we do not want to apply filter to the clinical mutations input files
  # run_maf_filter_unindexed:
  #   run: maf_filter.cwl
  #   in:
  #     sample: unindexed_samples
  #     maf_file:
  #       valueFrom: ${ return inputs.sample['maf_file']; }
  #     is_impact: is_impact
  #     argos_version_string: argos_version_string
  #   scatter: sample
  #   out: [ cbio_mutation_data_file ]

  # update the samples to use the new filtered maf files and output a single list of samples
  merge_samples_replace_mafs:
    in:
      samples: samples
        # source: [ samples, unindexed_samples ]
        # linkMerge: merge_flattened
      maf_files: run_maf_filter/cbio_mutation_data_file
        # source: [ run_maf_filter/cbio_mutation_data_file, run_maf_filter_unindexed/cbio_mutation_data_file ]
        # linkMerge: merge_flattened
    out: [ samples ]
    run:
      class: ExpressionTool
      inputs:
        samples:
          type:
            type: array
            items:
              type: record
              fields:
                maf_file: File
                sample_id: string
                normal_id: string
        maf_files: File[]
      outputs:
        samples:
          type:
            type: array
            items:
              type: record
              fields:
                maf_file: File
                sample_id: string
                normal_id: string
      # NOTE: in the line below `var i in inputs.samples`, `i` is an int representing the index position in the array `inputs.samples`
      # in Python it would look like ` x = ['a', 'b']; for i in range(len(x)): print(i, x[i]) `
      expression: "${
        var new_samples = [];

        for ( var i in inputs.samples ){
            new_samples.push({
              'sample_id': inputs.samples[i]['sample_id'],
              'normal_id': inputs.samples[i]['normal_id'],
              'maf_file': inputs.maf_files[i]
            });
          };

        return {'samples': new_samples};
        }"

  # run the fillout workflow
  run_samples_fillout:
    run: samples_fillout_workflow.cwl
    in:
      output_fname: fillout_output_fname
      exac_filter: exac_filter
      # samples: merge_samples_replace_mafs/samples
      samples:
        source: [ merge_samples_replace_mafs/samples, unindexed_samples ]
        linkMerge: merge_flattened
      bam_files:
        source: [ bam_files, run_indexer/bam_indexed ]
        linkMerge: merge_flattened
      ref_fasta: ref_fasta
    out: [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: run_samples_fillout/output_file
