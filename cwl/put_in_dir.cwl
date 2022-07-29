#!/usr/bin/env cwl-runner
cwlVersion: v1.0
doc: put all the items from a single list of Files or Directories into a new dir
class: ExpressionTool
id: put-in-dir

requirements:
  - class: InlineJavascriptRequirement

inputs:
  output_directory_name: string
  files:
    doc: a list of Files or Directories
    type:
      type: array
      items:
        - File
        - Directory
        - string
        - 'null'


outputs:

  directory:
    type: Directory

# This tool returns a Directory object,
# which holds all output files from the list
# of supplied input files
expression: |
  ${
    var output_files = [];
    var input_files = inputs.files.filter(single_file => String(single_file).toUpperCase() != 'NONE');


    for (var i = 0; i < inputs.files.length; i++) {
      if(input_files[i]){
        output_files.push(inputs.files[i]);
      }
    }

    return {
      'directory': {
        'class': 'Directory',
        'basename': inputs.output_directory_name,
        'listing': output_files
      }
    };
  }
