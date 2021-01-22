#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['reduce_sig_figs_seg.mean.py']
stdout: output.txt

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.01.1

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1

outputs:
  output_file:
    type: stdout
