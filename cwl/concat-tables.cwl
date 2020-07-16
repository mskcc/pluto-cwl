#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "concat-tables.py" ]
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:20.07.1
inputs:
  output_filename:
    type: string
    inputBinding:
      prefix: -o
      position: 1
  na_str:
    type: [ "null", string ]
    default: "NA"
    inputBinding:
      prefix: -n
      position: 2
  input_files:
    type: File[]
    inputBinding:
      position: 3

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
