#!/usr/bin/env cwl-runner
# bgzip -c "${vcf_sorted}" > "${vcf_gz}"
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bgzip" ]

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:igv-reports-1.0.1

stdout: $(inputs.output_filename)

inputs:
  input_file:
    type: File
    inputBinding:
      prefix: '-c'
      position: 1
  output_filename:
    type: string
    default: "output.gz"

outputs:
  output_file:
    type: stdout
