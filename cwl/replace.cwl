#!/usr/bin/env cwl-runner

# replace strings in file
# NOTE: in the future, do not write hard-coded CWL's like this, use a generic script or method

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ 'sed', 's/ILLOGICAL/NA/g' ]

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.3.4

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
