#!/usr/bin/env cwl-runner

# USAGE:
# $ filterUncalledMutations /path/to/data_mutations_extended.txt outputDir/

# NOTE:
# https://github.com/mskcc/pluto-cwl/issues/63
# https://github.com/mskcc/pluto-cwl/issues/84

cwlVersion: v1.2
class: CommandLineTool
baseCommand: ['bash', 'run.sh']
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix:21.5.2
  InitialWorkDirRequirement:
    listing:
      - entryname: run.sh
        entry: |-
          set -eu
          input_file="${ return inputs.input_file.path ; }"
          filterUncalledMutations "\${input_file}" .
          # output:
          # data_clinical_mutations.txt
          # data_mutations_uncalled.txt

inputs:
  input_file: File
  # output_filename:
  #   type: string
  #   default: "cases.txt"
outputs:
  called_file:
    type: File
    outputBinding:
      glob: data_clinical_mutations.txt
  uncalled_file:
    type: File
    outputBinding:
      glob: data_mutations_uncalled.txt
