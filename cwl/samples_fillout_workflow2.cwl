#!/usr/bin/env cwl-runner
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

steps:
  # convert all maf input files back to vcf
  maf2vcf:
    scatter: [ sample_id, maf_file ]
    scatterMethod: dotproduct
    in:
      sample_id: sample_ids
      maf_file: maf_files
      ref_fasta: ref_fasta
    out:
      [ output_file ]
    run:
      class: CommandLineTool
      baseCommand: ['bash', 'run.sh']
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:latest
        InitialWorkDirRequirement:
          listing:
            # NOTE: might need dos2unix for some that give errors ERROR: Your MAF uses CR line breaks, which we can't support. Please use LF or CRLF.
            # NOTE: might also need sanity check that maf has >1 line
            - entryname: run.sh
              entry: |-
                set -eu
                fasta="${ return inputs.ref_fasta.path; }"
                input_maf="${return inputs.maf_file.path;}"
                vcf="${ return inputs.sample_id + '.vcf'; }"
                vcf_sorted="${ return inputs.sample_id + '.sorted.vcf' }"
                vcf_sorted_gz="${ return inputs.sample_id + '.sorted.vcf.gz' }"
                # convert maf to vcf
                maf2vcf.pl --output-dir . --input-maf "\${input_maf}" --output-vcf "\${vcf}" --ref-fasta "\${fasta}"
                # sort the vcf
                bcftools sort --output-type v -o "\${vcf_sorted}" "\${vcf}"
                # archive the vcf
                bgzip -c "\${vcf_sorted}" > "\${vcf_sorted_gz}"
                # index the vcf
                tabix "\${vcf_sorted_gz}"
      # arguments:
      #   - valueFrom: $(inputs.maf_file)
      #     position: 1
      #     prefix: --input-maf
      #   - valueFrom: ${ return inputs.sample_id + '.vcf' }
      #     position: 2
      #     prefix: --output-vcf
      #   - valueFrom: $(inputs.ref_fasta)
      #     position: 3
      #     prefix: --ref-fasta
      inputs:
        sample_id: string
        maf_file: File
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
      outputs:
        output_file:
          type: File
          outputBinding:
            glob: ${ return inputs.sample_id + '.sorted.vcf.gz' }
          secondaryFiles:
            - .tbi

  # merge all the vcf files together to create a list of the target regions for all samples
  merge_vcfs:
    in:
      sample_ids: sample_ids
      vcf_gz_files: maf2vcf/output_file
    out:
      [ merged_vcf, merged_vcf_gz ]
    run:
      class: CommandLineTool
      baseCommand: ['bash', 'run.sh']
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:latest
        InitialWorkDirRequirement:
          listing:
            - entryname: run.sh
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
        sample_ids: string[]
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

  # run GetBaseCountsMultiSample on all the bam files against the target regions
  gbcms:
    in:
      sample_ids: sample_ids
      bam_files: bam_files
      targets_vcf: merge_vcfs/merged_vcf
      ref_fasta: ref_fasta
    out: [ output_file ]
    run:
      class: CommandLineTool
      baseCommand: [ "sh", "run.sh" ]
      requirements:
        InlineJavascriptRequirement: {}
        ResourceRequirement:
          coresMin: 4
        DockerRequirement:
          dockerPull: cmopipeline/getbasecountsmultisample:1.2.2
        InitialWorkDirRequirement:
          listing:
            - entryname: run.sh
              entry: |-
                set -eu
                # make the bam arg string;
                # --bam Sample5:Sample5.rg.md.abra.printreads.bam --bam Sample6:Sample6.rg.md.abra.printreads.bam
                bams_arg='${
                  var args1 = inputs.sample_ids.map((k, i) => [`${k}:${inputs.bam_files[i].path}`])
                  var args2 = args1.map( (a) => "--bam " + a )
                  return args2.join(" ") ;
                  }'
                fasta="${ return inputs.ref_fasta.path; }"
                vcf="${ return inputs.targets_vcf.path }"
                fillout_vcf="fillout.vcf"
                GetBaseCountsMultiSample --fasta "\${fasta}" --vcf "\${vcf}" --maq 20 --baq 20 --filter_improper_pair 0 --thread 4 --output "\${fillout_vcf}" \${bams_arg}
      inputs:
        sample_ids: string[]
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
        bam_files:
          type:
              type: array
              items: File
          secondaryFiles:
              - ^.bai
        targets_vcf:
          type: File
          # secondaryFiles:
          #   - .tbi
      outputs:
        output_file:
          type: File
          outputBinding:
            glob: fillout.vcf






outputs:
  gbcms_vcf:
    type: File
    outputSource: gbcms/output_file
  # sample_vcfs:
  #   type: File[]
  #   outputSource: maf2vcf/output_file
  # merged_vcf:
  #   type: File
  #   outputSource: merge_vcfs/output_file
