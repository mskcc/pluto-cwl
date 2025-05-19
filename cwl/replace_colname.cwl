#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["bash", "run.replace_colname.sh"]

requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: run.replace_colname.sh
        entry: |-
          set -eu
          old="$1"
          new="$2"
          input_file="$3"
          output_file="$4"
          awk -v old="\${old}" -v new="\${new}" 'NR==1 { gsub(old, new, $0); quit };1' "\${input_file}" > "\${output_file}"
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3

inputs:
  old_name:
    type: string
    inputBinding:
      position: 1
  new_name:
    type: string
    inputBinding:
      position: 2
  input_file:
    type: File
    inputBinding:
      position: 3
  output_filename:
    type: string
    default: output.tsv
    inputBinding:
      position: 4

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
