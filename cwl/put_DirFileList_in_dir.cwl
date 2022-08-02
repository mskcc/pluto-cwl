#!/usr/bin/env cwl-runner
cwlVersion: v1.0
doc: put all the items from a list of lists of Files or Directories into a new dir
class: ExpressionTool
id: put-DirFileList-dir

requirements:
  - class: InlineJavascriptRequirement

inputs:
  output_directory_name: string
  files:
    doc: a list of lists of Files or Directories
    type:
      type: array
      items:
        type: array
        items:
          - File
          - Directory
          - 'null'


outputs:
  directory:
    type: Directory

expression: |
  ${
    var output_files = [];

    for (var i in inputs.files) {
      if(String(inputs.files[i]).toUpperCase() != 'NONE'){
        for (var q in inputs.files[i]){
          if(String(inputs.files[i][q]).toUpperCase() != 'NONE'){
            output_files.push(inputs.files[i][q]);
          };
        };
      };
    };

    return {
      'directory': {
        'class': 'Directory',
        'basename': inputs.output_directory_name,
        'listing': output_files
      }
    };
  }
