#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bash", "run.mergebed.sh" ]

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.4.1
  InitialWorkDirRequirement:
    listing:
      - entryname: inputs_dir
        writable: true
        entry: "$({class: 'Directory', listing: inputs.bed_files})"
      - entryname: run.mergebed.sh
        # parallel --jobs 1 --xargs bedops -m {} ">" merged.{#}.bed ; <- chunks to largest size the system can handle but hit this error;
        # https://github.com/bedops/bedops/issues/249
        # so chunk it in groups of 1000 instead
        entry: |-
          set -euo pipefail
          find inputs_dir -type f | parallel --jobs 1 -n 1000 bedops -m {} ">" merged.{#}.bed
          find . -maxdepth 1 ! -path 'inputs_dir*' -type f -name "merged.*.bed" | parallel bedops -m {} ">" merged.bed

inputs:
  bed_files:
    type: File[]

outputs:
  output_file:
    type: File
    outputBinding:
      glob: merged.bed
