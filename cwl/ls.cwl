#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "ls", "-la" ]
doc: "CWL to save a copy of the execution Directory for debugging"

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.01.1

# this work but stages files in current dir
  # InitialWorkDirRequirement:
  #   listing: $(inputs.input_files)

# doesnt work...
  # InitialWorkDirRequirement:
  #   listing:
  #     - entry: "$({class: 'Directory', listing: inputs.input_files, basename: 'some_dir'})"
  #       entryname: some_dir

  InitialWorkDirRequirement:
    listing:
      - entryname: some_dir
        writable: true
        # entry: ${
        #   if ( inputs.input_files == null ) {
        #     return {class: 'Directory', listing: []};
        #   } else {
        #     return {class: 'Directory', listing: inputs.input_files};
        #   }
        # }
        entry: "$({class: 'Directory', listing: inputs.input_files})"


# this does not work
  # InitialWorkDirRequirement:
  #   listing:
  #     - class: Directory
  #       basename: some_dir
  #       listing: $(inputs.input_files)
      # - entryname: some_dir
      #   entry: $(inputs.input_files)




# doesnt work
  # InitialWorkDirRequirement:
  #   listing:
  #     - entry: $(inputs.input_files)
  #       # entryname: some_dir

# doesnt work
  # InitialWorkDirRequirement:
  #   listing:
  #     - class: Directory
  #       basename: some_dir
  #       listing: [$(inputs.input_files)]


    #
    # $(inputs.input_files)

# # doesnt work
#   InitialWorkDirRequirement:
#     listing:
#       - class: Directory
#         basename: some_dir
#         listing: []

      # - entryname: some_dir
      #   entry: $(inputs.input_files)
  #   listing:
  #     - entryname: run_facets_wrapper.sh
  #       entry: |-
  #         run-facets-wrapper.R --legacy-output TRUE --everything -D . --facets-lib-path /usr/local/lib/R/site-library $@ || touch failed.txt

stdout: ls.txt

inputs:
  input_files:
    type: File[]?
    # inputBinding:
    #   position: 4


outputs:
  output_file:
    type: stdout
