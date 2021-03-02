#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "maf2vcf.pl", "--output-dir", "." ]

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:igv-reports-1.0.1

inputs:
  maf_file:
    type: File
    inputBinding:
      prefix: '--input-maf'
      position: 1
  output_vcf_filename:
    type: string
    default: "output.vcf"
    inputBinding:
      prefix: '--output-vcf'
      position: 2
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
    inputBinding:
      prefix: '--ref-fasta'
      position: 3

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_vcf_filename)
