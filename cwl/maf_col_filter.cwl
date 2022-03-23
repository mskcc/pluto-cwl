#!/usr/bin/env cwl-runner
# CWL for script that removes extra columns from maf file

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['maf_col_filter.py']

# NOTE: bump this to the release version when its ready!
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.4.0

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
  output_filename:
    type: string
    inputBinding:
      position: 2

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
