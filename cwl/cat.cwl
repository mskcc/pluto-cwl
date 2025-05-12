#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ cat ]
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.4.1
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3
stdout: output.txt
inputs:
  input_files:
    type: File[]
    inputBinding:
      position: 1
outputs:
  output_file:
    type: stdout
