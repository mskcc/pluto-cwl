#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ cut ]
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.02.2
stdout: output.txt
inputs:
  field_indexes:
    type: string
    inputBinding:
      prefix: -f
      position: 1

  input_file:
    type: File
    inputBinding:
      position: 2
outputs:
  output_file:
    type: stdout
