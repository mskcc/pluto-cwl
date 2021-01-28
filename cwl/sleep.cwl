#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bash", "run.sh" ]

requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: run.sh
        entry: |-
          sleep 10

inputs:
  dummy_input:
    type: string

outputs: []
