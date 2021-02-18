#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "paste-col.py" ]
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.02.2
inputs:
  input_file:
    type: File
    inputBinding:
      prefix: -i
      position: 1
  output_filename:
    type: string
    inputBinding:
      prefix: -o
      position: 2
  header:
    type: string
    inputBinding:
      prefix: --header
      position: 3
  value:
    type: string
    inputBinding:
      prefix: --value
      position: 4

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
