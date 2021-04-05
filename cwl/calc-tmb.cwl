#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["calc-tmb.py", 'from-file']

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.03.1

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
  output_filename:
    type: string
    inputBinding:
      position: 2
  genome_coverage:
    type: string
    inputBinding:
      prefix: --genome-coverage
      position: 3
  normal_id:
    type: string
    inputBinding:
      prefix: --normal-id
      position: 4

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
