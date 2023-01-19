#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bash", "run.write_str.sh" ]

stdout: output.txt

requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: run.write_str.sh
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
