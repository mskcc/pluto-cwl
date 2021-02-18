#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
doc: "
Example CWL workflow that uses some advanced features
"

inputs:
  value:
    type: string
  sample_id:
    type: string

steps:
  # write the sample ID to a file
  make_file:
    run: write-str.cwl
    in:
      str: sample_id
    out:
      [ output_file ]

  # put a header line on the file
  add_header:
    run: add_header.cwl
    in:
      input_file: make_file/output_file
      header_str:
        valueFrom: ${ return "SampleID"; }
    out:
      [ output_file ]

  # add a second column to the file
  add_value_col:
    run: paste-col.cwl
    in:
      input_file: add_header/output_file
      output_filename:
        valueFrom: ${ return "output.tsv"; }
      header:
        valueFrom: ${ return "Value"; }
      value: value
    out:
      [ output_file ]

  # sleep for a few seconds so that we can see the system activity
  sleep:
    run: sleep.cwl
    in:
      dummy_input: value
    out:
      []

outputs:
  output_file:
    type: File
    outputSource: add_value_col/output_file
