#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bash", "sleep.sh" ]

requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: sleep.sh
        entry: |-
          sleep 5
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3

inputs:
  dummy_input:
    type: string

outputs: []
