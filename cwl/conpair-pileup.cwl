#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool
baseCommand:
  - python
  - /usr/bin/conpair/scripts/run_gatk_pileup_for_sample.py

id: conpair-pileup

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMin: 16000
    coresMin: 1
  DockerRequirement:
    dockerPull: mskcc/roslin-variant-conpair:0.3.3

doc: |
  None

inputs:

  ref:
    type: File
    inputBinding:
      prefix: --reference
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
      - .fai
      - ^.dict
  java_xmx:
    type:
    - 'null'
    - type: array
      items: string
    doc: set up java -Xmx parameter
    inputBinding:
      prefix: --xmx_java

  java_temp:
    type: ['null', string]
    doc: temporary directory to set -Djava.io.tmpdir
    inputBinding:
      prefix: --temp_dir_java

  gatk:
    type:
    - [File, string, "null"]
    inputBinding:
      prefix: --gatk

  markers_bed:
    type:
    - [File, string]
    inputBinding:
      prefix: --markers

  bam:
    type:
    - [File, string]
    inputBinding:
      prefix: --bam
    secondaryFiles:
      - ^.bai

  outfile:
    type:
    - string
    inputBinding:
      prefix: --outfile

outputs:
  out_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.outfile)
            return inputs.outfile;
          return null;
        }
