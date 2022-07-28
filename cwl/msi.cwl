#!/usr/bin/env cwl-runner

# example command;
#    msisensor msi \
#    -d msi_sites
#    -n normal_bam
#    -t tumor_bam
#    -tumor_sample_name tumor_sample_name
#    -normal_sample_name normal_sample_name
#    -o ${ return inputs.tumor_sample_name +"."+inputs.normal_sample_name+".msi.txt" }
#

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "msisensor", "msi" ]

requirements:
  DockerRequirement:
    dockerPull: mskcc/msisensor:0.2
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMin: 16000
    coresMin: 8

doc: Run msisensor on tumor-normal bams to differentiate MSI (microsatellite instable) samples from MSS (microsatellite stable) ones

# NOTE: next time use more verbose input labels, and only make input args for options we are using in the pipeline
inputs:
  threads:
    type: string
    default: "8"
    inputBinding:
      prefix: -b
  d:
    type:
    - string
    - File
    doc: homopolymer and microsatellites file
    inputBinding:
      prefix: -d

  n:
    type:
    - File
    doc: normal bam file
    secondaryFiles: ["^.bai"]
    inputBinding:
      prefix: -n

  t:
    type:
    - File
    doc: tumor bam file
    secondaryFiles: ["^.bai"]
    inputBinding:
      prefix: -t

  o:
    type: string
    doc: output distribution file
    inputBinding:
      prefix: -o

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.o)
  dis_file:
    type: File
    outputBinding:
      glob: $(inputs.o)_dis
  somatic_file:
    type: File
    outputBinding:
      glob: $(inputs.o)_somatic
