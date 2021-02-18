#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bash", "run.sh" ]

requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: run.sh
        # need a unique filename so that we can merge them later in a work dir without the names colliding. but also need parts of the name that are static so we can glob against it. but those static parts also need to be unique enough that we shouldnt have an input file matching that name
        # also we are just going to assume that input maf format has comments, a header starting with 'Hugo_Symbol', and ordered columns
        # NOTE: do not need to use the UUID here anymore since we can now rename files on staging
        entry: |-
          set -euo pipefail
          output_file=\$(python3 -c "import uuid; print('_maf2bed_merged.' + str(uuid.uuid4()) + '.bed')")
          grep -v '#' "$1" | grep -v 'Hugo' | cut -f5-7 | sort -V -k1,1 -k2,2n > "\$output_file"

inputs:
  maf_file:
    type: File
    inputBinding:
      position: 1

outputs:
  output_file:
    type: File
    outputBinding:
      glob: _maf2bed_merged.*.bed
