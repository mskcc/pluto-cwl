#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
doc: "
Version of concat-mafs.cwl that doesnt drop any columns or make columns blank
"

baseCommand: [
  "concat-tables.py",
  '--dir', # read input from directory of files
  '--na-str', '.', # using '' crashes GetBaseCountsMultiSample
  '--comments', # parse any header comments
  '--no-carriage-returns' # carriage returns output by default by Python csv.DictWriter will crash GetBaseCountsMultiSample
  ]

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.03.0
  InitialWorkDirRequirement:
    listing:
      - entryname: inputs_dir
        writable: true
        # put all the input files in a dir and rename them so the filenames do not collide
        entry: |-
            ${
              for (var i = 0; i < inputs.input_files.length; i++) {
                inputs.input_files[i].basename = inputs.input_files[i].basename + "." + i;
              }
              return {class: 'Directory', listing: inputs.input_files};
            }

arguments:
  - valueFrom: $(inputs.output_filename)
    position: 1
    prefix: -o
  - valueFrom: ${ return "inputs_dir" }
    position: 2

inputs:
  output_filename:
    type: string
    default: "output.maf"
  input_files:
    type: File[]

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
