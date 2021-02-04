#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "concat-tables.py", '--dir' ]
doc: "Concatenate all the table files provided, but put the input files in a dir first to try avoiding issues with gigantic CLI args when tons of files are passed"
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:latest
  InitialWorkDirRequirement:
    listing:
      - entryname: inputs_dir
        writable: true
        entry: "$({class: 'Directory', listing: inputs.input_files})"

arguments:
  - valueFrom: $(inputs.output_filename)
    position: 1
    prefix: -o
  - valueFrom: $(inputs.na_str)
    position: 2
    prefix: -n
  - valueFrom: $(inputs.comments)
    position: 3
    prefix: --comments
  - valueFrom: ${ return "inputs_dir" }
    position: 4

inputs:
  output_filename:
    type: string
  na_str:
    type: [ "null", string ]
    default: "NA"
  comments:
    type: [ "null", boolean ]
  input_files:
    type: File[]

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
