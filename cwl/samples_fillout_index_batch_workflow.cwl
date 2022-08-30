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
  sample_groups:
    doc: a nested list of lists of grouped samples
    type:
      type: array
      items:
        type: array
        items: "types.yml#FilloutIndexSample"
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
  # PREPROCESSING STEPS
  split_singleton_groups:
    doc: separate out any sample groups that had only 1 sample so they can be handled separately downstream
    in:
      sample_groups: sample_groups
    out: [ sample_groups, singleton_mafs ]
    run:
      class: ExpressionTool
      inputs:
        sample_groups:
          type:
            type: array
            items:
              type: array
              items: "types.yml#FilloutIndexSample"
      outputs:
        # singleton_samples: "types.yml#FilloutIndexSample[]"
        singleton_mafs: File[]
        sample_groups:
          type:
            type: array
            items:
              type: array
              items: "types.yml#FilloutIndexSample"
      expression: |
        ${
          var sample_groups = [];
          var singleton_mafs = [];

          for ( var i in inputs.sample_groups ){
            var group = inputs.sample_groups[i];

            if (group.length > 1) {
              sample_groups.push(group);
            } else if (group.length == 1){
              singleton_mafs.push(group[0].maf_file)
            }

          }

        return { 'sample_groups': sample_groups, "singleton_mafs": singleton_mafs };
        }

  # PRIMARY PROCESSING
  run_samples_fillout_index_workflow:
    doc: run the fillout workflow on each sample group with >1 sample
    scatter: samples
    # run: samples_fillout_index_workflow.cwl
    in:
      samples: split_singleton_groups/sample_groups
      ref_fasta: ref_fasta
      exac_filter: exac_filter
      is_impact: is_impact
      argos_version_string: argos_version_string
      fillout_output_fname: fillout_output_fname
    out: [ output_file, filtered_file, portal_file, uncalled_file ]
    run:
      class: Workflow
      id: run_samples_fillout_index_workflow
      label: run_samples_fillout_index_workflow
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
        is_impact:
          type: boolean
          default: True
        argos_version_string:
          type: [ "null", string ]
          default: "Unspecified"
        fillout_output_fname:
          type: string
          default: "fillout.maf"
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
      steps:
        fillout_index_prefilter:
          run: fillout_index_prefilter.cwl
          in:
            samples: samples
            is_impact: is_impact
            argos_version_string: argos_version_string
          out: [ samples ]
        run_samples_fillout:
          doc: run the fillout workflow on the samples
          run: samples_fillout_workflow.cwl
          in:
            output_fname: fillout_output_fname
            exac_filter: exac_filter
            samples: fillout_index_prefilter/samples
            ref_fasta: ref_fasta
          out: [ output_file, filtered_file, portal_file, uncalled_file ]





  # POST-PROCESSING - MERGE SAMPLE GROUPS FILES
  concat_output_files:
    run: concat-tables.cwl
    in:
      input_files: run_samples_fillout_index_workflow/output_file
      output_filename: fillout_output_fname
      na_str:
        valueFrom: ${ return ''; }
      comments:
        valueFrom: ${ return true; }
    out: [ output_file ]

  concat_filtered_file:
    run: concat-tables.cwl
    in:
      input_files: run_samples_fillout_index_workflow/filtered_file
      output_filename:
        valueFrom: ${ return 'output.filtered.maf'; }
      na_str:
        valueFrom: ${ return ''; }
      comments:
        valueFrom: ${ return true; }
    out: [ output_file ]

  concat_portal_file:
    run: concat-tables.cwl
    in:
      input_files: run_samples_fillout_index_workflow/portal_file
      output_filename:
        valueFrom: ${ return 'data_mutations_extended.txt'; }
      na_str:
        valueFrom: ${ return ''; }
      comments:
        valueFrom: ${ return true; }
    out: [ output_file ]

  concat_uncalled_file:
    run: concat-tables.cwl
    in:
      input_files: run_samples_fillout_index_workflow/uncalled_file
      output_filename:
        valueFrom: ${ return 'data_mutations_uncalled.txt'; }
      na_str:
        valueFrom: ${ return ''; }
      comments:
        valueFrom: ${ return true; }
    out: [ output_file ]

  # POST-PROCESSING - MERGE SINGLETON FILES
  concat_singletons:
    doc: merge the input maf files from each singleton sample
    in:
      singleton_mafs: split_singleton_groups/singleton_mafs
      concat_maf: concat_output_files/output_file
    out: [output_file]
    run:
      class: CommandLineTool
      baseCommand: [ "bash", "run.sh" ]
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:21.4.1
        InitialWorkDirRequirement:
          listing:
            - entryname: run.sh
              entry: |-
                set -eu
                # get a space-delim string of file paths
                singleton_mafs="${ return inputs.singleton_mafs.map((a) => a.path).join(' ') }"
                concat_maf="$(inputs.concat_maf.path)"
                output_maf="output.maf"
                # check how many singleton files were present
                num_singletons="\$(echo \${singleton_mafs} | wc -w )"
                # only try to concat if theres >0 singletons
                if [ "\${num_singletons}" -gt 0 ]; then
                concat-tables.py -o "\${output_maf}" --no-carriage-returns --comments --progress --na-str '' \${concat_maf} \${singleton_mafs}
                # otherwise just copy the input as the output
                else
                cp "\${concat_maf}" "\${output_maf}"
                fi

      inputs:
        singleton_mafs: File[]
        concat_maf: File
      outputs:
        output_file:
          type: File
          outputBinding:
            glob: output.maf



outputs:
  output_file:
    type: File
    outputSource: concat_singletons/output_file
  filtered_file:
    type: File
    outputSource: concat_filtered_file/output_file
  portal_file:
    type: File
    outputSource: concat_portal_file/output_file
  uncalled_file:
    type: File
    outputSource: concat_uncalled_file/output_file
