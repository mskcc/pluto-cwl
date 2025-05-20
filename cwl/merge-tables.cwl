#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "merge-tables.py" ]
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.4.1
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3

inputs:
  table1:
    type: File
    inputBinding:
      position: 1
  table2:
    type: File
    inputBinding:
      position: 2
  key1:
    type: string
    inputBinding:
      position: 3
      prefix: --key1
  key2:
    type: string
    inputBinding:
      position: 4
      prefix: --key2
  output_filename:
    type: string
    inputBinding:
      prefix: --output
      position: 5
  cBioPortal: # if the output table should have cBioPortal headers
    type: [ "null", boolean ]
    inputBinding:
      prefix: --cBioPortal
      position: 6

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
