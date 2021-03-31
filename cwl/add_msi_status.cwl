#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "add_msi_status.py" ]
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.03.0

inputs:
  input_filename:
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

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
