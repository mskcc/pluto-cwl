#!/usr/bin/env cwl-runner
# CWL for compiling a report from cBioPortal files
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "compile.R" ]
requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:reporting

arguments:
- valueFrom: "$(runtime.tmpdir)"
  position: 0
  prefix: '--intermediates'
  shellQuote: false

inputs:
  mutation_file:
    type: File
    inputBinding:
      position: 1
      prefix: '--mutations'
  samples_file:
    type: File
    inputBinding:
      position: 2
      prefix: '--samples'
  patients_file:
    type: File
    inputBinding:
      position: 3
      prefix: '--patients'
  output_filename:
    type: string
    default: report.html
    inputBinding:
      position: 4
      prefix: '--output_file'

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
