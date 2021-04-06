#!/usr/bin/env cwl-runner
# maf_filter.py "${maf}" "${config.version}" "${config.is_impact}" "${analysis_mut_file}" "${portal_file}"
# analysis_mut_file = "${config.project_id}.muts.maf"
#   portal_file = "data_mutations_extended.txt"

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['maf_filter.py']

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.3.2

inputs:
  maf_file:
    type: File
    inputBinding:
      position: 1
  argos_version_string:
    type: string
    inputBinding:
      prefix: '--version-string'
      position: 2
  is_impact:
    type: boolean
    inputBinding:
      prefix: '--is-impact'
      position: 3
  analysis_mutations_filename:
    type: ["null", string]
    default: "analysis.muts.maf"
    inputBinding:
      prefix: '--analyst-file'
      position: 4
  cbio_mutation_data_filename:
    type: ["null", string]
    default: "data_mutations_extended.txt"
    inputBinding:
      prefix: '--portal-file'
      position: 5
  rejected_file:
    type: [ "null", string ]
    default: "rejected.muts.maf"
    inputBinding:
      prefix: '--rejected-file'
      position: 6
  keep_rejects:
    type: boolean
    default: true
    inputBinding:
      prefix: '--keep-rejects'
      position: 7
outputs:
  analysis_mutations_file:
    type: File
    outputBinding:
      glob: $(inputs.analysis_mutations_filename)
  cbio_mutation_data_file:
    type: File
    outputBinding:
      glob: $(inputs.cbio_mutation_data_filename)
  rejected_file: # optional output
    type: File?
    outputBinding:
      glob: $(inputs.rejected_file)
