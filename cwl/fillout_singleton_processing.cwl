#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
id: fillout_singleton_processing
label: fillout_singleton_processing
doc: "
Handling for unmatched singleton sample files in the fillout workflow
"
requirements:
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
  - $import: types.yml

inputs:
  samples: "types.yml#FilloutIndexSample[]"
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
  exac_filter: # need this to resolve error in subworkflow: Anonymous file object must have 'contents' and 'basename' fields.
  # TODO: this needs the .tbi/.csi index file added!!
    type: File
    default:
      class: File
      path: /juno/work/ci/resources/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
  # these are needed for the filter script
  is_impact:
    type: boolean
    default: True
  argos_version_string:
    type: [ "null", string ]
    default: "Unspecified"

steps:
  fillout_index_prefilter:
    run: fillout_index_prefilter.cwl
    in:
      samples: samples
      is_impact: is_impact
      argos_version_string: argos_version_string
    out: [ samples ] # "types.yml#FilloutSample[]"

  fillout_pre_processing:
    run: fillout_pre_processing.cwl
    in:
      samples: fillout_index_prefilter/samples
      ref_fasta: ref_fasta
    out:
      [ sample_ids, clinical_sample_ids, bam_files, vcf_gz_files, merged_vcf, merged_vcf_gz ]




  # Need to make sure vcf has AD, FL_VF, and SRC fields for use in clinical filter !!
  # at this step, the vcf has these fields only;
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths of REF and ALT(s) in the order listed">
  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
  # https://samtools.github.io/hts-specs/VCFv4.3.pdf
  # https://samtools.github.io/bcftools/howtos/FAQ.html
  # https://github.com/samtools/bcftools/issues/810
  # https://samtools.github.io/bcftools/howtos/annotate.html
  fix_labels_and_merge_vcfs:
    doc: Need to update vcf to add AD, FL_VF, and SRC labels, to match downstream vcfs
    in:
      vcf: fillout_pre_processing/merged_vcf
    out: [ vcf ]
    run:
      class: CommandLineTool
      id: fix_labels_and_merge_vcfs
      label: fix_labels_and_merge_vcfs
      baseCommand: [ "bash", "run.sh" ]
      requirements:
        - class: DockerRequirement
          dockerPull: mskcc/helix_filters_01:21.4.1
        - class: InitialWorkDirRequirement
          listing:
            - entryname: run.sh
              entry: |-
                set -eux
                touch fix_labels_and_merge_vcfs
                # echo "$(inputs.vcf.path)"
                input_vcf="${ return inputs.vcf.path; }"
                output_vcf=updated.vcf
                echo "\${input_vcf}"
                # get the original header values, add some entries
                grep '##' "\${input_vcf}" > header
                echo '##INFO=<ID=SRC,Type=String,Number=.,Description="Source samples for the variant">' >> header
                echo '##FORMAT=<ID=FL_VF,Number=1,Type=Float,Description="Variant frequence (AD/DP)">' >> header

                # build the annotation table for SRC
                ANNOT_SRC=annot.SRC.txt
                printf '#CHROM\\tPOS\\tREF\\tALT\\tINFO\\n' > "\${ANNOT_SRC}"
                bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%SAMPLE,]\\n' -i 'AD>0'  "\${input_vcf}" >> "\${ANNOT_SRC}"
                bgzip -c "\${ANNOT_SRC}" > "\${ANNOT_SRC}".gz
                tabix -p vcf "\${ANNOT_SRC}".gz

                # creates the SRC column from the INFO field in annotations
                VCF_SRC=updated.SRC.vcf
                bcftools annotate \\
                --header-lines header \\
                --annotations "\${ANNOT_SRC}".gz \\
                -c CHROM,POS,REF,ALT,SRC \\
                --output "\${VCF_SRC}" \\
                --output-type v \\
                "\${input_vcf}"

                # create the annotation table for FL_VF, use '.' as empty value
                ANNOT_FLVF=annot.FLVF.txt
                printf '#CHROM\\tPOS\\tREF\\tALT\\tFL_VF\\n' > "\${ANNOT_FLVF}"
                bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t.\n'  "\${input_vcf}" >> "\${ANNOT_FLVF}"
                bgzip -c "\${ANNOT_FLVF}" > "\${ANNOT_FLVF}".gz
                tabix -p vcf "\${ANNOT_FLVF}".gz

                # make samples list
                bcftools query -l "\${input_vcf}" > samples.txt

                # add annotation for FL_VF
                VCF_FLVF=updated.FLVF.vcf
                bcftools annotate \\
                --header-lines header \\
                --annotations "\${ANNOT_FLVF}".gz \\
                -c CHROM,POS,REF,ALT,FORMAT/FL_VF \\
                -S <(head -1 samples.txt) \\
                --output "\${VCF_FLVF}" \\
                --output-type v \\
                "\${VCF_SRC}"
                cp "\${VCF_SRC}" "\${output_vcf}"
      inputs:
        vcf: File
      outputs:
        vcf:
          type: File
          outputBinding:
            glob: updated.vcf



  fillout_post_processing:
    run: fillout_post_processing.cwl
    in:
      samples: fillout_index_prefilter/samples
      fillout_vcf: fix_labels_and_merge_vcfs/vcf
      clinical_sample_ids: fillout_pre_processing/clinical_sample_ids
      ref_fasta: ref_fasta
      exac_filter: exac_filter
    out: [ output_file, filtered_file, portal_file, uncalled_file ]

outputs:
  output_file:
    type: File
    outputSource: fillout_post_processing/output_file
  filtered_file:
    type: File
    outputSource: fillout_post_processing/filtered_file
  portal_file:
    type: File
    outputSource: fillout_post_processing/portal_file
  uncalled_file:
    type: File
    outputSource: fillout_post_processing/uncalled_file
