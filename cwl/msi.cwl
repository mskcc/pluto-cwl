#!/usr/bin/env cwl-runner
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
    doc: threads number for parallel computing
    type: string
    default: "8"
    # NOTE: no performance gains seen when using >8 threads
    inputBinding:
      prefix: -b

  microsatellites_file:
    type:
    - string
    - File
    doc: homopolymer and microsatellites file
    inputBinding:
      prefix: -d

  normal_bam:
    type:
    - File
    doc: normal bam file
    secondaryFiles: ["^.bai"]
    inputBinding:
      prefix: -n

  tumor_bam:
    type:
    - File
    doc: tumor bam file
    secondaryFiles: ["^.bai"]
    inputBinding:
      prefix: -t

  output_filename:
    type: string
    doc: output distribution file
    default: msi.txt
    inputBinding:
      prefix: -o

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
  dis_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)_dis
  somatic_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)_somatic


# example command;
#    msisensor msi \
#    -d msi_sites
#    -n normal_bam
#    -t tumor_bam
#    -tumor_sample_name tumor_sample_name
#    -normal_sample_name normal_sample_name
#    -o ${ return inputs.tumor_sample_name +"."+inputs.normal_sample_name+".msi.txt" }
#

#
# Program: msisensor (homopolymer and miscrosatelite analysis using bam files)
# Version: v0.2
# Author: Beifang Niu && Kai Ye
#
# Usage:   msisensor <command> [options]
#
# Key commands:
#
#  scan            scan homopolymers and miscrosatelites
#  msi             msi scoring
#
#
#
# Singularity> msisensor msi
#
# Usage:  msisensor msi [options]
#
#        -d   <string>   homopolymer and microsates file
#        -n   <string>   normal bam file
#        -t   <string>   tumor  bam file
#        -o   <string>   output distribution file
#
#        -e   <string>   bed file, optional
#        -f   <double>   FDR threshold for somatic sites detection, default=0.05
#        -c   <int>      coverage threshold for msi analysis, WXS: 20; WGS: 15, default=20
#        -r   <string>   choose one region, format: 1:10000000-20000000
#        -l   <int>      mininal homopolymer size, default=5
#        -p   <int>      mininal homopolymer size for distribution analysis, default=10
#        -m   <int>      maximal homopolymer size for distribution analysis, default=50
#        -q   <int>      mininal microsates size, default=3
#        -s   <int>      mininal microsates size for distribution analysis, default=5
#        -w   <int>      maximal microstaes size for distribution analysis, default=40
#        -u   <int>      span size around window for extracting reads, default=500
#        -b   <int>      threads number for parallel computing, default=1
#        -x   <int>      output homopolymer only, 0: no; 1: yes, default=0
#        -y   <int>      output microsatellite only, 0: no; 1: yes, default=0
#
#        -h   help
