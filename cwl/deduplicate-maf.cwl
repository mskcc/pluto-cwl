#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bash", "run.sh" ]
doc: "Remove duplicate lines from the maf file

Expected columns are in this order;

1  Hugo_Symbol
2  Center
3  NCBI_Build
4  Chromosome
5  Start_Position
6  End_Position
7  Variant_Classification
8  Reference_Allele
9  Tumor_Seq_Allele1
10  Tumor_Seq_Allele2
11  n_alt_count
12  Matched_Norm_Sample_Barcode
13  t_alt_count
14  t_ref_count
15  n_ref_count
16  Tumor_Sample_Barcode

NOTE: actually this order might be wrong but all that matters is that Chromosome and Start_Position are columns 4 and 5
"

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.3.2
  InitialWorkDirRequirement:
    listing:
      - entryname: run.sh
        # need to get comments
        # then get header
        # then sort the maf rows
        entry: |-
          set -euo pipefail
          grep '#' "$1" > comments || touch comments
          grep '^Hugo_Symbol' "$1" > header || touch header
          grep -v '^Hugo_Symbol' "$1" | grep -v '#' | sort -u > body.dedup
          sort -V -k4,4 -k5,5n body.dedup > body.sorted
          cat comments > dedup.maf
          cat header >> dedup.maf
          cat body.sorted >> dedup.maf
inputs:
  input_file:
    type: File
    inputBinding:
      position: 1

outputs:
  output_file:
    type: File
    outputBinding:
      glob: dedup.maf
