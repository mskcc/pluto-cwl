#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "ls", "-la", "some_dir" ]
doc: "CWL to save a copy of the execution Directory for debugging"

requirements:
  InlineJavascriptRequirement: {}
  # DockerRequirement:
  #   dockerPull: mskcc/helix_filters_01:21.03.0
  InitialWorkDirRequirement:
    listing:
      - entryname: some_dir
        writable: true
        entry: |-
            ${
              for (var i = 0; i < inputs.input_files.length; i++) {
                inputs.input_files[i].basename = inputs.input_files[i].basename + "." + i;
              }
              return {class: 'Directory', listing: inputs.input_files};
            }

stdout: ls.txt

inputs:
  input_files:
    type: File[]

outputs:
  output_file:
    type: stdout
