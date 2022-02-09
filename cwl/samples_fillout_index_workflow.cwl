#!/usr/bin/env cwl-runner

cwlVersion: v1.1
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
          sample_id: string # must match sample ID used inside maf file
          normal_id: string
          prefilter: boolean # should the sample input maf file be filtered before including for fillout
          sample_type: string # "research" or "clinical" to dictate downstream handling and germline filtering
          maf_file: File
          bam_file: File

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
  # bam file indexing subworkflow
  index_all_bams:
    in:
      samples: samples
    out: [ samples ]
    run:
      class: Workflow
      inputs:
        samples:
          type:
            type: array
            items:
              type: record
              fields:
                sample_id: string
                normal_id: string
                prefilter: boolean
                sample_type: string
                maf_file: File
                bam_file: File
      outputs:
        samples:
          outputSource: update_sample_bams/samples
          type:
            type: array
            items:
              type: record
              fields:
                maf_file: File
                sample_id: string
                normal_id: string
                prefilter: boolean
                sample_type: string
                bam_file:
                  type: File
                  secondaryFiles:
                      - ^.bai
      steps:
        # get a list of just the bam files for downstream processing
        create_bam_list:
          in:
            samples: samples
          out:
            [ bam_files ]
          run:
            class: ExpressionTool
            inputs:
              samples:
                type:
                  type: array
                  items:
                    type: record
                    fields:
                      bam_file: File
            outputs:
              bam_files:
                type:
                    type: array
                    items: File
                secondaryFiles:
                    - ^.bai
            expression: "${
              var bam_files = [];
              for ( var i in inputs.samples ){
                  bam_files.push(inputs.samples[i]['bam_file']);
                };
              return {'bam_files': bam_files};
              }"
        # index each bam file
        run_indexer:
          run: index_bam.cwl
          in:
            bam: create_bam_list/bam_files
          scatter: bam
          out: [ bam_indexed ]
        # update each sample's bam file with the new indexed file from the list
        update_sample_bams:
          in:
            samples: samples
            bam_files: run_indexer/bam_indexed
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
                      sample_id: string
                      normal_id: string
                      prefilter: boolean
                      sample_type: string
                      maf_file: File
                      bam_file: File
              bam_files: File[]
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
                      prefilter: boolean
                      sample_type: string
                      bam_file:
                        type: File
                        secondaryFiles:
                            - ^.bai
            # NOTE: in the line below `var i in inputs.samples`, `i` is an int representing the index position in the array `inputs.samples`
            # in Python it would look like ` x = ['a', 'b']; for i in range(len(x)): print(i, x[i]) `
            expression: "${
              var new_samples = [];

              for ( var i in inputs.samples ){
                  new_samples.push({
                    'sample_id': inputs.samples[i]['sample_id'],
                    'normal_id': inputs.samples[i]['normal_id'],
                    'prefilter': inputs.samples[i]['prefilter'],
                    'sample_type': inputs.samples[i]['sample_type'],
                    'maf_file': inputs.samples[i]['maf_file'],
                    'bam_file': inputs.bam_files[i]
                  });
                };
              console.log(new_samples);
              return {'samples': new_samples};
              }"





  # some samples need their input maf files filtered ahead of time; apply the maf_filter cBioPortal filters
  apply_prefilter:
    in:
      samples: index_all_bams/samples
      is_impact: is_impact
      argos_version_string: argos_version_string
    out: [ samples ]
    run:
      class: Workflow
      inputs:
        is_impact:
          type: boolean
          default: True
        argos_version_string:
          type: [ "null", string ]
          default: "Unspecified"
        samples:
          type:
            type: array
            items:
              type: record
              fields:
                sample_id: string
                normal_id: string
                prefilter: boolean
                sample_type: string
                maf_file: File
                bam_file:
                  type: File
                  secondaryFiles:
                      - ^.bai
      outputs:
        samples:
          outputSource: run_maf_filter_on_samples/samples
          type:
            type: array
            items:
              type: record
              fields:
                maf_file: File
                sample_id: string
                normal_id: string
                prefilter: boolean
                sample_type: string
                bam_file:
                  type: File
                  secondaryFiles:
                      - ^.bai
      steps:
        # need to separate the samples that require prefilter from the ones that do not
        split_sample_groups:
          in:
            samples: samples
          out: [ samples_need_filter, samples_no_filter ]
          run:
            class: ExpressionTool
            inputs:
              samples:
                type:
                  type: array
                  items:
                    type: record
                    fields:
                      sample_id: string
                      normal_id: string
                      prefilter: boolean
                      sample_type: string
                      maf_file: File
                      bam_file:
                        type: File
                        secondaryFiles:
                            - ^.bai
            outputs:
              samples_need_filter:
                type:
                  type: array
                  items:
                    type: record
                    fields:
                      sample_id: string
                      normal_id: string
                      prefilter: boolean
                      sample_type: string
                      maf_file: File
                      bam_file:
                        type: File
                        secondaryFiles:
                            - ^.bai
              samples_no_filter:
                type:
                  type: array
                  items:
                    type: record
                    fields:
                      sample_id: string
                      normal_id: string
                      prefilter: boolean
                      sample_type: string
                      maf_file: File
                      bam_file:
                        type: File
                        secondaryFiles:
                            - ^.bai
            # also consider: String(x).toLowerCase() == "true"
            expression: "${
              var samples_no_filter = [];
              var samples_need_filter = [];

              for ( var i in inputs.samples ){
                if ( inputs.samples[i]['prefilter'] === true ) {
                    samples_need_filter.push(inputs.samples[i]);
                  } else {
                    samples_no_filter.push(inputs.samples[i]);
                  }
              };

              console.log(samples_no_filter);
              console.log(samples_need_filter);

              return {
                  'samples_need_filter': samples_need_filter,
                  'samples_no_filter': samples_no_filter
                };
              }"

        # run the maf_filter on the samples_need_filter and return the cBioPortal files for each
        # need a subworkflow for this because maf_filter.cwl takes a single input file
        # then need to update the samples channel with the new output files
        run_maf_filter_on_samples:
          in:
            samples: split_sample_groups/samples_need_filter
            samples_no_filter: split_sample_groups/samples_no_filter
            is_impact: is_impact
            argos_version_string: argos_version_string
          out: [ samples ]
          run:
            class: Workflow
            inputs:
              is_impact:
                type: boolean
                default: True
              argos_version_string:
                type: [ "null", string ]
                default: "Unspecified"
              samples:
                type:
                  type: array
                  items:
                    type: record
                    fields:
                      sample_id: string
                      normal_id: string
                      prefilter: boolean
                      sample_type: string
                      maf_file: File
                      bam_file:
                        type: File
                        secondaryFiles:
                            - ^.bai
              samples_no_filter:
                type:
                  type: array
                  items:
                    type: record
                    fields:
                      sample_id: string
                      normal_id: string
                      prefilter: boolean
                      sample_type: string
                      maf_file: File
                      bam_file:
                        type: File
                        secondaryFiles:
                            - ^.bai
            outputs:
              samples:
                outputSource: update_sample_mafs/samples
                type:
                  type: array
                  items:
                    type: record
                    fields:
                      maf_file: File
                      sample_id: string
                      normal_id: string
                      prefilter: boolean
                      sample_type: string
                      bam_file:
                        type: File
                        secondaryFiles:
                            - ^.bai
            steps:
              # need to separate out a channel of just maf_files in an array
              get_maf_files:
                in:
                  samples: samples
                out: [ maf_files ]
                run:
                  class: ExpressionTool
                  inputs:
                    samples:
                      type:
                        type: array
                        items:
                          type: record
                          fields:
                            sample_id: string
                            normal_id: string
                            prefilter: boolean
                            sample_type: string
                            maf_file: File
                            bam_file:
                              type: File
                              secondaryFiles:
                                  - ^.bai
                  outputs:
                    maf_files: File[]
                  expression: "${
                    var maf_files = [];
                    for ( var i in inputs.samples ){
                      maf_files.push( inputs.samples[i]['maf_file'] );
                    }
                    console.log(maf_files);
                    return {'maf_files': maf_files};
                    }"
              run_maf_filter:
                run: maf_filter.cwl
                in:
                  maf_file: get_maf_files/maf_files
                  is_impact: is_impact
                  argos_version_string: argos_version_string
                scatter: maf_file
                out: [ cbio_mutation_data_file ]
              # update the maf file for each sample that was filtered; merge in the unfiltered samples as well
              update_sample_mafs:
                in:
                  samples: samples
                  samples_no_filter: samples_no_filter
                  maf_files: run_maf_filter/cbio_mutation_data_file
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
                            sample_id: string
                            normal_id: string
                            prefilter: boolean
                            sample_type: string
                            maf_file: File
                            bam_file:
                              type: File
                              secondaryFiles:
                                  - ^.bai
                    maf_files: File[]
                    samples_no_filter:
                      type:
                        type: array
                        items:
                          type: record
                          fields:
                            sample_id: string
                            normal_id: string
                            prefilter: boolean
                            sample_type: string
                            maf_file: File
                            bam_file:
                              type: File
                              secondaryFiles:
                                  - ^.bai
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
                            prefilter: boolean
                            sample_type: string
                            bam_file:
                              type: File
                              secondaryFiles:
                                  - ^.bai
                  expression: "${
                    var new_samples = [];
                    for ( var i in inputs.samples ){
                        new_samples.push({
                          'sample_id': inputs.samples[i]['sample_id'],
                          'normal_id': inputs.samples[i]['normal_id'],
                          'prefilter': inputs.samples[i]['prefilter'],
                          'sample_type': inputs.samples[i]['sample_type'],
                          'bam_file': inputs.samples[i]['bam_file'],
                          'maf_file': inputs.maf_files[i]
                        });
                      };
                      for ( var i in inputs.samples_no_filter ){
                          new_samples.push({
                            'sample_id': inputs.samples_no_filter[i]['sample_id'],
                            'normal_id': inputs.samples_no_filter[i]['normal_id'],
                            'prefilter': inputs.samples_no_filter[i]['prefilter'],
                            'sample_type': inputs.samples_no_filter[i]['sample_type'],
                            'bam_file': inputs.samples_no_filter[i]['bam_file'],
                            'maf_file': inputs.samples_no_filter[i]['maf_file']
                          });
                        };
                    console.log(new_samples);
                    return {'samples': new_samples};
                    }"


  # run the fillout workflow
  run_samples_fillout:
    run: samples_fillout_workflow.cwl
    in:
      output_fname: fillout_output_fname
      exac_filter: exac_filter
      samples: apply_prefilter/samples
      ref_fasta: ref_fasta
    out: [ output_file ]


outputs:
  output_file:
    type: File
    outputSource: run_samples_fillout/output_file
