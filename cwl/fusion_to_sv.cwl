#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["fusion_to_sv_converter.py"]
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.6.0
inputs:
  fusion_file:
    type: File
    inputBinding:
      prefix: --fusion_file
      position: 1
  output_filename:
    type: string
    inputBinding:
      prefix: --sv_file
      position: 2

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
