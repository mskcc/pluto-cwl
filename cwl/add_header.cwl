#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["bash", "run.add_header.sh"]
stdout: output.txt
requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: run.add_header.sh
        entry: |-
          echo "$1"
          cat "$2"
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3
inputs:
  header_str:
    type: string
    inputBinding:
      position: 1
  input_file:
    type: File
    inputBinding:
      position: 2

outputs:
  output_file:
    type: stdout
