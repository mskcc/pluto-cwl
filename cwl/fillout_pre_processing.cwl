#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
id: fillout_pre_processing
label: fillout_pre_processing
doc: "

"

requirements:
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
  - $import: types.yml

inputs:
  samples: "types.yml#FilloutSample[]" # must be a list of samples; list must have at  least one entry 
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
  create_samples_list:
    id: create_samples_list
    doc: creates a list of sample_ids out of the samples record array for downstream use
    in:
      samples: samples
    out: [ sample_ids, num_samples ]
    run:
      class: ExpressionTool
      id: create_samples_list
      label: create_samples_list
      inputs:
        samples: "types.yml#FilloutSample[]"
      outputs:
        sample_ids: string[]
        num_samples: int
      # NOTE: in the line below `var i in inputs.samples`, `i` is an int representing the index position in the array `inputs.samples`
      # in Python it would look like ` x = ['a', 'b']; for i in range(len(x)): print(i, x[i]) `
      expression: |
        ${
        var sample_ids = [];
        var num_samples = 0; // need to output the number of samples because I dont know how to check the list length in the CWL when statement
        for ( var i in inputs.samples ){
            sample_ids.push(inputs.samples[i]['sample_id']);
          };
        return {'sample_ids': sample_ids, 'num_samples': sample_ids.length};
        }

  create_clinical_samples_list:
    doc: creates a list of clinical sample_ids out of the samples record array for downstream use
    in:
      samples: samples
    out: [ clinical_sample_ids ]
    run:
      class: ExpressionTool
      id: create_clinical_samples_list
      label: create_clinical_samples_list
      inputs:
        samples: "types.yml#FilloutSample[]"
      outputs:
        clinical_sample_ids: string[]
      expression: |
        ${
          var sample_ids = [];
          for ( var i in inputs.samples ){
              if (inputs.samples[i]['sample_type'] === 'clinical' ){
                sample_ids.push(inputs.samples[i]['sample_id']);
              };
            };
          return {'clinical_sample_ids': sample_ids};
        }

  create_bam_list:
    doc: gets a list of just the bam files for downstream process
    in:
      samples: samples
    out:
      [ bam_files ]
    run:
      class: ExpressionTool
      id: create_bam_list
      label: create_bam_list
      inputs:
        samples: "types.yml#FilloutSample[]"
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


  maf2vcf:
    run: fillout_maf2vcf.cwl
    doc: converts all maf input files back to vcf for downstream processing
    scatter: sample
    in:
      sample: samples
      sample_id:
        valueFrom: ${ return inputs.sample['sample_id']; }
      maf_file:
        valueFrom: ${ return inputs.sample['maf_file']; }
      ref_fasta: ref_fasta
    out:
      [ output_file, output_vcf ] # .sorted.vcf.gz + .tbi, .sorted.vcf

  merge_vcfs:
    doc: merge all the vcf files together to create the list of target regions for GetBaseCountsMultiSample for fillout
    in:
      num_samples: create_samples_list/num_samples
      sample_ids: create_samples_list/sample_ids
      vcf_gz_files: maf2vcf/output_file
    out:
      [ merged_vcf, merged_vcf_gz ]
    when: $( inputs.num_samples > 1 )
    run:
      class: CommandLineTool
      id: merge_vcfs
      label: merge_vcfs
      baseCommand: ['bash', 'run.merge_vcfs.sh']
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:21.4.1
        InitialWorkDirRequirement:
          listing:
            - entryname: run.merge_vcfs.sh
              entry: |-
                set -eu
                # get a comma-delim string of the sample names
                samples_arg="${ return inputs.sample_ids.join(',') }"
                # get a space-delim string of vcf file paths
                vcf_files="${ return inputs.vcf_gz_files.map((a) => a.path).join(' ') }"
                merged_vcf="merged.vcf"
                merged_vcf_gz="merged.vcf.gz"
                bcftools merge --force-samples --output-type v \${vcf_files} | \\
                bcftools sort | \\
                bcftools view -s "\${samples_arg}" > "\${merged_vcf}"
                bgzip -c "\${merged_vcf}" > "\${merged_vcf_gz}"
                tabix "\${merged_vcf_gz}"
      inputs:
        sample_ids: 
          type: string[]
          doc: sample id's MUST match what are used in the header columns in the .vcf files!
        vcf_gz_files:
          type:
              type: array
              items: File
          secondaryFiles:
              - .tbi
      outputs:
        merged_vcf:
          type: File
          outputBinding:
            glob: merged.vcf
        merged_vcf_gz:
          type: File
          outputBinding:
            glob: merged.vcf.gz
          secondaryFiles:
            - .tbi

  pass_singleton_vcf:
    doc: this passes along the original sample .vcf if only one sample was provided and no merge_vcf would have been produced
    in:
      num_samples: create_samples_list/num_samples
      vcf_gz_files: maf2vcf/output_file
      vcf_files: maf2vcf/output_vcf
    when: $( inputs.num_samples == 1 )
    out: [ output_vcf_gz, output_vcf ]
    run:
      class: ExpressionTool
      id: pass_singleton_vcf
      label: pass_singleton_vcf
      inputs:
        vcf_gz_files: File[]
        vcf_files: File[]
      outputs:
        output_vcf_gz:
          type: File
          secondaryFiles:
            - .tbi
        output_vcf: File
      expression: |
        ${
          return {
            'output_vcf_gz': inputs.vcf_gz_files[0],
            'output_vcf': inputs.vcf_files[0]
            };
        }




outputs: 
  sample_ids:
    type: string[]
    outputSource: create_samples_list/sample_ids
  clinical_sample_ids:
    type: string[]
    outputSource: create_clinical_samples_list/clinical_sample_ids
  bam_files:
    type: File[]
    outputSource: create_bam_list/bam_files
  vcf_gz_files:
    type: File[]
    outputSource: maf2vcf/output_file
  merged_vcf:
    type: File
    outputSource: [ merge_vcfs/merged_vcf, pass_singleton_vcf/output_vcf ]
    pickValue: first_non_null
  merged_vcf_gz:
    type: File
    outputSource: [ merge_vcfs/merged_vcf_gz, pass_singleton_vcf/output_vcf_gz ]
    pickValue: first_non_null
