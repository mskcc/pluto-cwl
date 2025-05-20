#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "ls", "-la" ]
doc: "CWL to save a copy of the execution Directory for debugging"

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3

stdout: ls.txt

inputs: []

outputs:
  output_file:
    type: stdout
