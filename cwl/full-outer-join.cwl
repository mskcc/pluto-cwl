#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "full-outer-join.R" ]
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:dev
inputs:
  table1:
    type: File
    inputBinding:
      position: 1
  table2:
    type: File
    inputBinding:
      position: 2
  join_key:
    type: string
    inputBinding:
      position: 3
  output_filename:
    type: string
    inputBinding:
      position: 4

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
