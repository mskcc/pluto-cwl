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
  sites_maf: File # single input maf file that will be used as the regions to use in IGV
  maf_files: File[] # list of all maf files to list as tracks in IGV
  bam_files:
    type:
      type: array
      items: File
    secondaryFiles: [^.bai]
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
  # need to convert the sites maf to vcf format since that is how we have the igv-report.cwl configured; this could be .bed in the future
  convert_sites_to_vcf:
    run: maf2vcf_gz_workflow.cwl
    in:
      maf_file: sites_maf
      ref_fasta: ref_fasta
    out:
      [ output_file ]

  convert_mafs_to_vcf:
    run: maf2vcf_gz_workflow.cwl
    scatter: maf
    in:
      maf: maf_files
      maf_file:
        valueFrom: ${ return inputs.maf; }
      ref_fasta: ref_fasta
    out:
      [ output_file ] # [ "output.vcf.gz", ... ]

  make_report:
    run: igv-reports.cwl
    in:
      sites: convert_sites_to_vcf/output_file
      vcf_gz_files: convert_mafs_to_vcf/output_file
      bam_files: bam_files
      ref_fasta: ref_fasta
    out:
      [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: make_report/output_file
