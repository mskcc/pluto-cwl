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
    out: [ sample_groups, singleton_mafs, singleton_samples, num_singletons ]
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
        singleton_samples: "types.yml#FilloutIndexSample[]"
        singleton_mafs: File[]
        num_singletons: int
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
          var singleton_samples = [];
          var num_singletons = 0; // need to output the number of singletons because I dont know how to check the list length in the CWL when statement

          for ( var i in inputs.sample_groups ){
            var group = inputs.sample_groups[i];

            if (group.length > 1) {
              sample_groups.push(group);
            } else if (group.length == 1){
              singleton_samples.push(group[0]);
              singleton_mafs.push(group[0].maf_file);
            };
          };

        var res = { 'sample_groups': sample_groups, 
          "singleton_mafs": singleton_mafs, 
          "singleton_samples": singleton_samples,
          "num_singletons": singleton_samples.length };

        return res;
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




  # # POST-PROCESSING - MERGE SINGLETON FILES
  # # NOTE: need to get FL_VF into the singleton .vcf !!!
  singleton_processing:
    run: fillout_singleton_processing.cwl
    in:
      num_singletons: split_singleton_groups/num_singletons
      samples: split_singleton_groups/singleton_samples
      ref_fasta: ref_fasta
      exac_filter: exac_filter
      is_impact: is_impact
      argos_version_string: argos_version_string
    when: $( inputs.num_singletons > 0 )
    # all outputs will be null if this step is skipped
    out: [ output_file, filtered_file, portal_file, uncalled_file ]
  
  concat_singleton_files:
    in:
      output_mafs: 
        source: [ concat_output_files/output_file, singleton_processing/output_file ] 
        pickValue: all_non_null
      filtered_mafs: 
        source: [ concat_filtered_file/output_file, singleton_processing/filtered_file ] 
        pickValue: all_non_null
      portal_mafs: 
        source: [ concat_portal_file/output_file, singleton_processing/portal_file ] 
        pickValue: all_non_null
      uncalled_mafs: 
        source: [ concat_uncalled_file/output_file, singleton_processing/uncalled_file ] 
        pickValue: all_non_null
    out: [ output_maf, filtered_maf, portal_maf, uncalled_maf ]
    run:
      class: CommandLineTool
      id: concat_singleton_files
      label: concat_singleton_files
      baseCommand: ["bash", "run.concat_singleton_files.sh"]
      requirements: 
        - class: DockerRequirement
          dockerPull: mskcc/helix_filters_01:21.4.1
        - class: InitialWorkDirRequirement
          listing:
            - entryname: run.concat_singleton_files.sh
              entry: |-
                set -eu
                num_output_mafs="${ return inputs.output_mafs.length; }"
                num_filtered_mafs="${ return inputs.filtered_mafs.length; }"
                num_portal_mafs="${ return inputs.portal_mafs.length; }"
                num_uncalled_mafs="${ return inputs.uncalled_mafs.length; }"
                if [ \${num_output_mafs} -gt 0 ]; then 
                  # # get a space-delim string of file paths
                  output_mafs="${ return inputs.output_mafs.map((a) => a.path).join(' ') }"
                  concat-tables.py -o "output.maf" --no-carriage-returns --comments --na-str '' \${output_mafs}
                else 
                  # idk what to do here really ???
                  touch output.maf
                fi

                if [ \${num_filtered_mafs} -gt 0 ]; then 
                  filtered_mafs="${ return inputs.filtered_mafs.map((a) => a.path).join(' ') }"
                  concat-tables.py -o "output.filtered.maf" --no-carriage-returns --comments --na-str '' \${filtered_mafs}
                else 
                  touch output.filtered.maf
                fi

                if [ \${num_portal_mafs} -gt 0 ]; then 
                  portal_mafs="${ return inputs.portal_mafs.map((a) => a.path).join(' ') }"
                  concat-tables.py -o "data_mutations_extended.txt" --no-carriage-returns --comments --na-str '' \${portal_mafs}
                else 
                  touch data_mutations_extended.txt
                fi

                if [ \${num_uncalled_mafs} -gt 0 ]; then 
                  uncalled_mafs="${ return inputs.uncalled_mafs.map((a) => a.path).join(' ') }"
                  concat-tables.py -o "data_mutations_uncalled.txt" --no-carriage-returns --comments --na-str '' \${uncalled_mafs}
                else 
                  touch data_mutations_uncalled.txt
                fi
      inputs:
        output_mafs: File[]
        filtered_mafs: File[]
        portal_mafs: File[]
        uncalled_mafs: File[]
      outputs:
        output_maf:
          type: File
          outputBinding:
            glob: output.maf
        filtered_maf:
          type: File
          outputBinding:
            glob: output.filtered.maf
        portal_maf:
          type: File
          outputBinding:
            glob: data_mutations_extended.txt
        uncalled_maf:
          type: File
          outputBinding:
            glob: data_mutations_uncalled.txt



outputs:
  output_file:
    type: File
    outputSource: concat_singleton_files/output_maf
  filtered_file:
    type: File
    outputSource: concat_singleton_files/filtered_maf
  portal_file:
    type: File
    outputSource: concat_singleton_files/portal_maf
  uncalled_file:
    type: File
    outputSource: concat_singleton_files/uncalled_maf
