#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: fillout_maf2vcf
label: fillout_maf2vcf

doc: converts all maf input files back to vcf for downstream processing
# NOTE: This is important; do NOT try to do complex manipulations on maf format file, do it on vcf format instead

baseCommand: ['bash', 'fillout_maf2vcf.run.sh']
requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: mskcc/helix_filters_01:21.4.1
  - class: InitialWorkDirRequirement
    listing:
    # NOTE: might need dos2unix for some that give errors ERROR: Your MAF uses CR line breaks, which we can't support. Please use LF or CRLF.
    # NOTE: might also need sanity check that maf has >1 line
    - entryname: fillout_maf2vcf.run.sh
      entry: |-
        set -eu
        fasta="${ return inputs.ref_fasta.path; }"
        input_maf="${return inputs.maf_file.path;}"
        vcf="${ return inputs.sample_id + '.vcf'; }"
        vcf_sorted="${ return inputs.sample_id + '.sorted.vcf' }"
        vcf_sorted_gz="${ return inputs.sample_id + '.sorted.vcf.gz' }"
        sample_id="${ return inputs.sample_id }"
        # convert maf to vcf
        num_lines="\$(wc -l < \${input_maf})"
        more_than_five="\$(( \${num_lines} > 5))"
        if [ "\${more_than_five}" == "1" ]
        then
          maf2vcf.pl --output-dir . --input-maf "\${input_maf}" --output-vcf "\${vcf}" --ref-fasta "\${fasta}"
        else
          cat << EOF > "\${vcf}"
        ##fileformat=VCFv4.2
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##FORMAT=<ID=FL_DP,Number=1,Type=Integer,Description="Total depth">
        ##FORMAT=<ID=FL_RD,Number=1,Type=Integer,Description="Depth matching reference (REF) allele">
        ##FORMAT=<ID=FL_AD,Number=1,Type=Integer,Description="Depth matching alternate (ALT) allele">
        ##FORMAT=<ID=FL_VF,Number=1,Type=Float,Description="Variant frequence (AD/DP)">
        ##FORMAT=<ID=FL_DPP,Number=1,Type=Integer,Description="Depth on postitive strand">
        ##FORMAT=<ID=FL_DPN,Number=1,Type=Integer,Description="Depth on negative strand">
        ##FORMAT=<ID=FL_RDP,Number=1,Type=Integer,Description="Reference depth on postitive strand">
        ##FORMAT=<ID=FL_RDN,Number=1,Type=Integer,Description="Reference depth on negative strand">
        ##FORMAT=<ID=FL_ADP,Number=1,Type=Integer,Description="Alternate depth on postitive strand">
        ##FORMAT=<ID=FL_ADN,Number=1,Type=Integer,Description="Alternate depth on negative strand">
        ##FORMAT=<ID=DPF,Number=1,Type=Integer,Description="Total fragment depth">
        ##FORMAT=<ID=RDF,Number=1,Type=Float,Description="Fragment depth matching reference (REF) allele">
        ##FORMAT=<ID=ADF,Number=1,Type=Float,Description="Fragment depth matching alternate (ALT) allele">
        ##contig=<ID=1>
        ##contig=<ID=4>
        ##contig=<ID=5>
        ##contig=<ID=10>
        ##contig=<ID=16>
        ##contig=<ID=19>
        ##contig=<ID=X>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths of REF and ALT(s) in the order listed">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
        ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
        ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
        ##bcftools_annotateVersion=1.9+htslib-1.9
        ##bcftools_annotateCommand=annotate --header-lines merged.header.txt --annotations /var/lib/cwl/stgcdd0cd9e-22d6-4b54-adf1-3504bbb3198f/merged.vcf.gz --columns FORMAT/GT,FORMAT/AD,FORMAT/DP,INFO/AC,INFO/AN --output fillout.merged.vcf --output-type v fillout.relabel.vcf.gz; Date=Wed Mar 29 17:10:46 2023
        ##INFO=<ID=SRC,Type=String,Number=.,Description="Source samples for the variant">
        ##bcftools_annotateCommand=annotate --header-lines fillout.merged.header.txt --annotations fillout.merged.annotation.txt.gz -c CHROM,POS,REF,ALT,SRC --output fillout.merged.sources.vcf --output-type v fillout.merged.vcf; Date=Wed Mar 29 17:10:47 2023
        ##bcftools_filterVersion=1.9+htslib-1.9
        ##bcftools_filterCommand=filter -e 'AD[@clinical_samples.txt:*]='.' && FL_VF[@clinical_samples.txt]>0.1' fillout.merged.sources.vcf; Date=Wed Mar 29 17:10:57 2023
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	\${sample_id}
        EOF
        fi
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
    doc: .vcf.gz file with .tbi index
    type: File
    outputBinding:
      glob: ${ return inputs.sample_id + '.sorted.vcf.gz' }
    secondaryFiles:
      - .tbi
  output_vcf:
    type: File
    outputBinding:
      glob: ${ return inputs.sample_id + '.sorted.vcf' }
