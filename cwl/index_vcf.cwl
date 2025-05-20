#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "tabix", "-p", "vcf" ]

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:igv-reports-1.0.1
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3
  # tabix does not work unless you stage the input file here
  InitialWorkDirRequirement:
    listing: [ $(inputs.input_file) ]

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1

outputs:
  indexed_file:
    type: File
    outputBinding:
      glob: $(inputs.input_file.basename)
    secondaryFiles: [".tbi"]
