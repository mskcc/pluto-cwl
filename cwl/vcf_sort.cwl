#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bcftools", "sort" ]

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:igv-reports-1.0.1
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3

inputs:
  vcf_file:
    type: File
    inputBinding:
      position: 2
  output_filename:
    type: string
    default: "output.sorted.vcf"
    inputBinding:
      prefix: '--output-file'
      position: 1

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
