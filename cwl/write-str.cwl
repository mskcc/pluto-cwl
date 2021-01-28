#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bash", "run.sh" ]

stdout: output.txt

requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: run.sh
        entry: |-
          echo "$1"

inputs:
  str:
    type: string
    inputBinding:
      position: 1

outputs:
  output_file:
    type: stdout
