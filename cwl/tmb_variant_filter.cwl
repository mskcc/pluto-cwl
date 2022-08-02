#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["tmb_variant_filter.py"]

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.4.1

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
  output_filename:
    type: string
    default: tmb_variants.maf
    inputBinding:
      position: 2

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
