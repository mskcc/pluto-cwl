#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "env" ]
doc: "CWL to save a copy of the execution environment for debugging"

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.03.1

stdout: env.container.txt

inputs: []

outputs:
  output_file:
    type: stdout
