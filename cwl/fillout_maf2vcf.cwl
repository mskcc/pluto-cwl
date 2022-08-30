#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: fillout_maf2vcf
label: fillout_maf2vcf

doc: converts all maf input files back to vcf for downstream processing
# NOTE: This is important; do NOT try to do complex manipulations on maf format file, do it on vcf format instead

baseCommand: ['bash', 'run.sh']
requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: mskcc/helix_filters_01:21.4.1
  - class: InitialWorkDirRequirement
    listing:
    # NOTE: might need dos2unix for some that give errors ERROR: Your MAF uses CR line breaks, which we can't support. Please use LF or CRLF.
    # NOTE: might also need sanity check that maf has >1 line
    - entryname: run.sh
      entry: |-
        set -eu
        fasta="${ return inputs.ref_fasta.path; }"
        input_maf="${return inputs.maf_file.path;}"
        vcf="${ return inputs.sample_id + '.vcf'; }"
        vcf_sorted="${ return inputs.sample_id + '.sorted.vcf' }"
        vcf_sorted_gz="${ return inputs.sample_id + '.sorted.vcf.gz' }"
        # convert maf to vcf
        maf2vcf.pl --output-dir . --input-maf "\${input_maf}" --output-vcf "\${vcf}" --ref-fasta "\${fasta}"
        # sort the vcf
        bcftools sort --output-type v -o "\${vcf_sorted}" "\${vcf}"
        # archive the vcf
        bgzip -c "\${vcf_sorted}" > "\${vcf_sorted_gz}"
        # index the vcf
        tabix "\${vcf_sorted_gz}"
inputs:
  sample_id: string
  maf_file: File
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
outputs:
  output_file:
    type: File
    outputBinding:
      glob: ${ return inputs.sample_id + '.sorted.vcf.gz' }
    secondaryFiles:
      - .tbi