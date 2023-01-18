#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
id: fillout_index_prefilter
label: fillout_index_prefilter
doc: "
Apply .bam indexing
and optional .maf pre-filtering
for fillout samples
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
    type: "types.yml#FilloutMafOptionalNoIndexSample[]" # "types.yml#FilloutIndexSample[]"
  # these are needed for the filter script
  is_impact:
    type: boolean
    default: True

  argos_version_string:
    type: [ "null", string ]
    default: "Unspecified"

steps:
  index_bam:
    doc: create a .bai index file for all incoming .bam files
    scatter: sample
    in:
      sample: samples
    out: [ sample ]
    run:
      class: CommandLineTool
      baseCommand: [ "bash", "run.index_bam.sh" ]
      requirements:
        ResourceRequirement:
          coresMin: 4
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:samtools-1.9
        InitialWorkDirRequirement:
          listing:
            - $(inputs.sample['bam_file']) # NOTE: I think this is causing all input bam files to get copied when using Toil??
            - entryname: run.index_bam.sh
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
        sample: "types.yml#FilloutMafOptionalNoIndexSample" # "types.yml#FilloutIndexSample"
      outputs:
        sample:
          type: "types.yml#FilloutMafOptionalIndexedSample"
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
          type: "types.yml#FilloutMafOptionalIndexedSample[]"
      outputs:
        samples_need_filter:
          type: "types.yml#FilloutIndexedSample[]"
        samples_no_filter:
          type: "types.yml#FilloutMafOptionalIndexedSample[]"
      # also consider: String(x).toLowerCase() == "true"
      expression: |
        ${
        var samples_no_filter = [];
        var samples_need_filter = [];

        for ( var i in inputs.samples ){
          // if it has a maf file and prefilter is true
          if ( inputs.samples[i]["prefilter"] === true && inputs.samples[i]["maf_file"] != null) {
              samples_need_filter.push(inputs.samples[i]);
            } else {
              samples_no_filter.push(inputs.samples[i]);
            }
        };

        return {
            'samples_need_filter': samples_need_filter,
            'samples_no_filter': samples_no_filter
          };
        }

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


  # convert FilloutIndexedSample to FilloutSample by removing some extraneous fields
  convert_sample_types:
    doc: Convert the CWL sample objects from FilloutIndexedSample type to FilloutSample for downstream processes
    in:
      # samples:
      #   source: [ split_sample_groups/samples_no_filter, prefilter_workflow/sample ]
      #   linkMerge: merge_flattened
      samples_no_filter: split_sample_groups/samples_no_filter
      samples_filtered: prefilter_workflow/sample
    out: [ samples ]
    run:
      class: ExpressionTool
      inputs:
        samples_no_filter: "types.yml#FilloutMafOptionalIndexedSample[]"
        samples_filtered: "types.yml#FilloutIndexedSample[]"
      outputs:
        samples:
          type: "types.yml#FilloutMafOptionalSample[]" # "types.yml#FilloutSample[]"
      expression: |
        ${
          var new_samples = [];

          for ( var i in inputs.samples_no_filter ){
            var d = inputs.samples_no_filter[i];
            delete d['prefilter'];
            new_samples.push(d);
          }

          for ( var i in inputs.samples_filtered ){
            var d = inputs.samples_filtered[i];
            delete d['prefilter'];
            new_samples.push(d);
          }

        return { 'samples': new_samples };
        }

outputs:
  samples:
    type: "types.yml#FilloutMafOptionalSample[]" # "types.yml#FilloutSample[]"
    outputSource: convert_sample_types/samples