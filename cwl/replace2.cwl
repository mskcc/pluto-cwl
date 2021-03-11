#!/usr/bin/env cwl-runner

# replace strings in file

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ 'sed', 's/%/MSI_SCORE/g' ]

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.02.2

stdout: $(inputs.output_filename)

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
  output_filename:
    type: ["null", string]
    default: "output.txt"

outputs:
  output_file:
    type: stdout
