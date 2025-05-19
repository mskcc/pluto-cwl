#!/usr/bin/env cwl-runner

# USAGE:
# $ filterUncalledMutations /path/to/data_mutations_extended.txt outputDir/

# NOTE:
# https://github.com/mskcc/pluto-cwl/issues/63
# https://github.com/mskcc/pluto-cwl/issues/84

cwlVersion: v1.2
class: CommandLineTool
baseCommand: ['bash', 'run.filterUncalledMutations.sh']
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix:21.5.3
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3
  InitialWorkDirRequirement:
    listing:
      - entryname: run.filterUncalledMutations.sh
        entry: |-
          set -eu
          input_file="${ return inputs.input_file.path ; }"
          filterUncalledMutations "\${input_file}" --output-dir .
          # output:
          # data_mutations_extended.txt
          # data_mutations_uncalled.txt

inputs:
  input_file: File
outputs:
  called_file:
    type: File
    outputBinding:
      glob: data_mutations_extended.txt
  uncalled_file:
    type: File
    outputBinding:
      glob: data_mutations_uncalled.txt
