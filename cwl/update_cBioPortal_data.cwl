#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "update_cBioPortal_data.py" ]

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:20.07.3

inputs:
  # required for every invocation
  subcommand:
    type: string
    inputBinding:
      position: 1
  input_file:
    type: File
    inputBinding:
      prefix: '--input'
      position: 2
  output_filename:
    type: string
    inputBinding:
      prefix: '--output'
      position: 3
  facets_txt:
    type: File
    inputBinding:
      prefix: '--facets-txt'
      position: 4

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
