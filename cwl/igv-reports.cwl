#!/usr/bin/env cwl-runner

# create_report \
# --standalone \
# --flanking 100 \
# --info-columns GT AD DP \
# --tracks "${vcf_gz}" \
# Sample5.rg.md.abra.printreads.bam \
# Sample6.rg.md.abra.printreads.bam \
# --output igv."${vcf_sorted}".html \
# "${vcf_gz}" \
# b37.fasta

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "create_report" ]

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:igv-reports-1.0.1
  ResourceRequirement:
    ramMin: 8000
    coresMin: 3

arguments:
  - valueFrom: ${ return inputs.vcf_gz_files.map((a) => a.path).join(' ') + ' ' + inputs.bam_files.map((a) => a.path).join(' '); }
    position: 1
    prefix: --tracks
    shellQuote: false
  - valueFrom: ${ return "100"; }
    position: 2
    prefix: --flanking
    shellQuote: false
  - valueFrom: ${ return "GT AD DP"; }
    position: 3
    prefix: --info-columns
    shellQuote: false
  - valueFrom: $(inputs.output_filename)
    prefix: --output
  - valueFrom: ${ return "--standalone"; }
    position: 4
    shellQuote: false
  - valueFrom: ${ return inputs.sites.path + ' ' + inputs.ref_fasta.path; }
    position: 5
    shellQuote: false

inputs:
  sites: # a .vcf.gz file with .tbi index ; TODO: this could also be a bed file I think, need a way to specify that
    type: File
    secondaryFiles: [.tbi]
  vcf_gz_files:
    type:
      type: array
      items: File
    secondaryFiles: [.tbi]
  bam_files:
    type:
      type: array
      items: File
    secondaryFiles: [^.bai]
  ref_fasta:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
      - .fai
      - ^.dict
  output_filename:
    type: string
    default: "igv.html"

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
