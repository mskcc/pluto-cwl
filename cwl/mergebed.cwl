#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bash", "run.sh" ]

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:latest
  InitialWorkDirRequirement:
    listing:
      - entryname: inputs_dir
        writable: true
        entry: "$({class: 'Directory', listing: inputs.bed_files})"
      - entryname: run.sh
        entry: |-
          set -euo pipefail
          find inputs_dir -type f | parallel --jobs 1 --xargs bedops -m {} ">" merged.{#}.bed
          find . -maxdepth 1 ! -path 'inputs_dir*' -type f -name "merged.*.bed" | parallel bedops -m {} ">" merged.bed

inputs:
  bed_files:
    type: File[]

outputs:
  output_file:
    type: File
    outputBinding:
      glob: merged.bed
