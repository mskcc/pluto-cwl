#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
baseCommand: ['bash', 'run.sh']
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix:21.4.2
  InitialWorkDirRequirement:
    listing:
      - entryname: run.sh
        entry: |-
          set -eu
          # get a comma-delim string of the sample names
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
