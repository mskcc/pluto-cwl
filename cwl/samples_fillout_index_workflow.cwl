#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
doc: "
Wrapper to run bam indexing on all bams before submitting for samples fillout
Also includes steps to pre-filter some maf input files
"
requirements:
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
  - $import: types.yml

inputs:
  samples:
    type: "types.yml#FilloutIndexSample[]"
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

  index_bam:
    doc: create a .bai index file for all incoming .bam files
    scatter: sample
    in:
      sample: samples
    out: [ sample ]
    run:
      class: CommandLineTool
      baseCommand: [ "bash", "run.sh" ]
      requirements:
        ResourceRequirement:
          coresMin: 4
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:samtools-1.9
        InitialWorkDirRequirement:
          listing:
            - $(inputs.sample['bam_file'])
            - entryname: run.sh
              entry: |-
                set -eu
                # sample.bam
                input_bam="$(inputs.sample['bam_file'].basename)"
                # sample.bam.bai
                default_bai="\${input_bam}.bai"
                # sample.bai
                extra_bai="\${input_bam%.*}.bai"
                samtools index -@ 4 "\${input_bam}"
                cp "\${default_bai}" "\${extra_bai}"

      inputs:
        sample: "types.yml#FilloutIndexSample"
      outputs:
        sample:
          type: "types.yml#FilloutIndexedSample"
          outputBinding:
            outputEval: ${
              var ret = inputs.sample;
              ret['bam_file']['secondaryFiles'] = [{"class":"File", "path":runtime.outdir + "/" + inputs.sample["bam_file"].nameroot + ".bai"}];
              return ret;
              }

  split_sample_groups:
    doc: separate out the samples that need to be prefiltered in downstream steps
    in:
      samples: index_bam/sample
    out: [ samples_need_filter, samples_no_filter ]
    run:
      class: ExpressionTool
      inputs:
        samples:
          type: "types.yml#FilloutIndexedSample[]"
      outputs:
        samples_need_filter:
          type: "types.yml#FilloutIndexedSample[]"
        samples_no_filter:
          type: "types.yml#FilloutIndexedSample[]"
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

        return {
            'samples_need_filter': samples_need_filter,
            'samples_no_filter': samples_no_filter
          };
        }"


  prefilter_workflow:
    doc: apply the prefiltering steps to applicable samples (usually Argos sample mafs since clinical mafs are often filtered already)
    scatter: sample
    in:
      sample: split_sample_groups/samples_need_filter
      is_impact: is_impact
      argos_version_string: argos_version_string
    out: [ sample ]
    run:
      class: Workflow
      inputs:
        is_impact:
          type: boolean
          default: True
        argos_version_string:
          type: [ "null", string ]
          default: "Unspecified"
        sample:
          type: "types.yml#FilloutIndexedSample"
      outputs:
        sample:
          outputSource: update_sample_mafs/sample
          type: "types.yml#FilloutIndexedSample"
      steps:
        run_maf_filter:
          doc: run the maf filer script on the sample maf file
          run: maf_filter.cwl
          in:
            sample: sample
            maf_file:
              valueFrom: $(inputs.sample['maf_file'])
            is_impact: is_impact
            argos_version_string: argos_version_string
          out: [ cbio_mutation_data_file ]
        update_sample_mafs:
          doc: update the sample entry to use the new filtered maf instead of the original one
          in:
            sample: sample
            maf_file: run_maf_filter/cbio_mutation_data_file
          out: [ sample ]
          run:
            class: ExpressionTool
            inputs:
              sample:
                type: "types.yml#FilloutIndexedSample"
              maf_file: File
            outputs:
              sample:
                type: "types.yml#FilloutIndexedSample"
            expression: |
              ${
              var new_sample = inputs.sample ;
              new_sample['maf_file'] = inputs.maf_file ;
              return { 'sample': new_sample };
              }



  convert_sample_types:
    doc: Convert the CWL sample objects from FilloutIndexedSample type to FilloutSample for downstream processes
    in:
      samples:
        source: [ split_sample_groups/samples_no_filter, prefilter_workflow/sample ]
        linkMerge: merge_flattened
    out: [ samples ]
    run:
      class: ExpressionTool
      inputs:
        samples:
          type: "types.yml#FilloutIndexedSample[]"
      outputs:
        samples:
          type: "types.yml#FilloutSample[]"
      expression: |
        ${
          var new_samples = [];
          for ( var i in inputs.samples ){
            var d = inputs.samples[i];
            delete d['prefilter'];
            new_samples.push(d);
          }
        return { 'samples': new_samples };
        }

  run_samples_fillout:
    doc: run the fillout workflow on the samples
    run: samples_fillout_workflow.cwl
    in:
      output_fname: fillout_output_fname
      exac_filter: exac_filter
      samples: convert_sample_types/samples
      ref_fasta: ref_fasta
    out: [ output_file, filtered_file, portal_file, uncalled_file ]


outputs:
  output_file:
    type: File
    outputSource: run_samples_fillout/output_file
  filtered_file:
    type: File
    outputSource: run_samples_fillout/filtered_file
  portal_file:
    type: File
    outputSource: run_samples_fillout/portal_file
  uncalled_file:
    type: File
    outputSource: run_samples_fillout/uncalled_file
