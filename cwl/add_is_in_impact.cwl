#!/usr/bin/env cwl-runner
# add_is_in_impact.py --input_file abcd.maf --output_file abcd_is_in_IMPACT_col_added.maf --IMPACT_file IMPACT505_b37_targets.bed

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['add_is_in_impact.py', '--include-assay']

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.03.1

inputs:
  input_file:
    type: File
    inputBinding:
      prefix: '--input_file'
      position: 1
  output_filename:
    type: string
    inputBinding:
      prefix: '--output_file'
      position: 2
  IMPACT_file: 
    type: File
    inputBinding:
      prefix: '--IMPACT_file'
      position: 3

outputs:
  IMPACT_col_added_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
