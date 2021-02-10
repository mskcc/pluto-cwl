#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
doc: "
Get the base counts in a .bam file for all the position in the supplied maf file
"
baseCommand: [
  "GetBaseCountsMultiSample",
  "--omaf",
  "--maq", "20",
  "--baq", "20",
  "--filter_improper_pair", "0"
  ]
#   --fasta b37.fasta \
# --maf Proj_08390_G.cols_subset.muts.maf \
# --output output.Proj_08390_G.cols_subset.muts.maf \
# --bam Sample6:Sample6.rg.md.abra.printreads.bam \
# --bam Sample5:Sample5.rg.md.abra.printreads.bam
# --thread

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: cmopipeline/getbasecountsmultisample:1.2.2

arguments:
  - valueFrom: $(inputs.output_filename)
    position: 1
    prefix: --output
  - valueFrom: $(inputs.maf)
    position: 2
    prefix: --maf
  - valueFrom: $(inputs.ref_fasta)
    position: 3
    prefix: --fasta
  - valueFrom: ${ return inputs.sample_id + ':' + inputs.bam.path; }
    position: 4
    prefix: --bam

inputs:
  sample_id:
    type: string
  ref_fasta:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
      - .fai
      - ^.dict
  bam:
    type: File
    secondaryFiles:
        - ^.bai
    doc: "Bam file for the sample that we want to check the base counts for"
  maf:
    type: File
    doc: "Maf file with positions that we want to check in the bam file"
  output_filename:
    type: string
    default: "fillout.maf"

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
