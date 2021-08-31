#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "full-outer-join.R" ]
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:R-3.5.1
inputs:
  table1:
    type: File
    inputBinding:
      position: 1
  table2:
    type:
      - 'null'
      - File[] 
    inputBinding:
      position: 2
      prefix: '--t2'
  join_key:
    type: string
    inputBinding:
      position: 3
      prefix: '-k'
  output_filename:
    type: string
    default: output.tsv
    inputBinding:
      position: 4
      prefix: '-o'

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
