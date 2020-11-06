#!/usr/bin/env cwl-runner
# CWL for script to add the variant allele frequency column for tumor sample to the maf file

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['add_af.py']

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:20.11.1

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
  output_filename:
    type: string
    inputBinding:
      position: 2

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
