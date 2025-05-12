#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
baseCommand: ['bash', 'run.updateCaseList.sh']
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix:21.5.1
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3
  InitialWorkDirRequirement:
    listing:
      - entryname: run.updateCaseList.sh
        entry: |-
          set -eu
          # get a comma-delim string of the sample names ; returns an empty string if there are no sample id's ; updateCaseList ignores empty strings
          samples_arg="${ return inputs.sample_ids.join(',') ; }"
          input_file="${ return inputs.case_list.path ; }"
          output_file="${ return inputs.output_filename }"

          updateCaseList "\${input_file}" "\${samples_arg}" > "\${output_file}"

inputs:
  sample_ids: string[]
  case_list: File
  output_filename:
    type: string
    default: "cases.txt"
outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
