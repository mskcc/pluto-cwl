#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
id: samples_fillout_workflow
label: samples_fillout_workflow
doc: "
Workflow to run GetBaseCountsMultiSample fillout on a number of samples, each with their own bam and maf files
"

requirements:
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
  - $import: types.yml

inputs:
  samples: "types.yml#FilloutSample[]" # type: array
  output_fname:
    type: [ 'null', string ]
    default: "output.maf"
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
  exac_filter:
    type: File
    secondaryFiles:
      - .tbi
    default:
      class: File
      path: /juno/work/ci/resources/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
  # NOTE: VEP cache dir is too large to use as an input item, need to pack it inside the container instead ?
  # vep_cache:
  #   type: Directory
  #   default:
  #     class: Directory
  #     path: /juno/work/ci/resources/vep/cache

steps:
  # PRE-PROCESSING
  fillout_pre_processing:
    run: fillout_pre_processing.cwl
    in:
      samples: samples
      ref_fasta: ref_fasta
    out:
      [ sample_ids, clinical_sample_ids, bam_files, vcf_gz_files, merged_vcf, merged_vcf_gz ]
  

  # PRIMARY FILLOUT PROCESSING STARTS HERE
  # TODO: convert this to a `scatter` step that runs per-sample in parallel, then merge the outputs otherwise we will hit the command line arg length issues eventually
  # NOTE: maybe do not do this ^^^ because currently its useful to have a multi-sample vcf output, need to investigate this more
  gbcms:
    doc: run GetBaseCountsMultiSample on all the bam files against the target regions (the merged vcf from all samples)
    in:
      sample_ids: fillout_pre_processing/sample_ids
      bam_files: fillout_pre_processing/bam_files
      targets_vcf: fillout_pre_processing/merged_vcf
      ref_fasta: ref_fasta
    out: [ output_file ]
    run:
      class: CommandLineTool
      id: gbcms
      label: gbcms
      baseCommand: [ "sh", "run.sh" ]
      requirements:
        InlineJavascriptRequirement: {}
        ResourceRequirement:
          coresMin: 8
        DockerRequirement:
          dockerPull: cmopipeline/getbasecountsmultisample:1.2.2
        InitialWorkDirRequirement:
          listing:
            - entryname: run.sh
              entry: |-
                set -eu
                # make the bam arg string; looks like this:
                # --bam Sample5:Sample5.rg.md.abra.printreads.bam --bam Sample6:Sample6.rg.md.abra.printreads.bam
                bams_arg='${
                  var args1 = inputs.sample_ids.map((k, i) => [`${k}:${inputs.bam_files[i].path}`])
                  var args2 = args1.map( (a) => "--bam " + a )
                  return args2.join(" ") ;
                  }'
                fasta="${ return inputs.ref_fasta.path; }"
                vcf="${ return inputs.targets_vcf.path }"
                fillout_vcf="fillout.vcf"
                GetBaseCountsMultiSample --fasta "\${fasta}" --vcf "\${vcf}" --maq 20 --baq 20 --filter_improper_pair 0 --thread 8 --output "\${fillout_vcf}" \${bams_arg}
      inputs:
        sample_ids: string[]
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
        bam_files:
          type:
              type: array
              items: File
          secondaryFiles:
              - ^.bai
        targets_vcf:
          type: File
      outputs:
        output_file:
          type: File
          outputBinding:
            glob: fillout.vcf

  fix_labels_and_merge_vcfs:
    doc: clean up, re-annotate, and re-header the fillout vcf, also merge against the fillout targets vcf in order to add a SRC sample source tag to tell which samples each variant originated in, and pull back in each variants original quality metrics
    in:
      fillout_vcf: gbcms/output_file
      merged_vcf: fillout_pre_processing/merged_vcf
      merged_vcf_gz: fillout_pre_processing/merged_vcf_gz
      ref_fasta: ref_fasta
    out: [ fillout_sources_vcf ]
    run:
      class: CommandLineTool
      id: fix_labels_and_merge_vcfs
      label: fix_labels_and_merge_vcfs
      baseCommand: [ "bash", "run.sh" ]
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:21.4.1
        InitialWorkDirRequirement:
          listing:
            - entryname: run.sh
              entry: |-
                set -eu
                # need to rename some conflicting sample INFO tags; prepend FL for "fillout"
                # change these tags: DP:RD:AD:VF:DPP:DPN:RDP:RDN:ADP:ADN
                # pre-pend them all with 'FL_' for fillout so they do not collide with downstream labels
                fillout_vcf="${ return inputs.fillout_vcf.path; }"
                fillout_relabel_vcf="fillout.relabel.vcf"
                fillout_relabel_vcf_gz="fillout.relabel.vcf.gz"
                sed -e 's|\\([=[:space:]]\\)DP\\([,:[:space:]]\\)|\\1FL_DP\\2|' \\
                -e 's|\\([=:[:space:]]\\)RD\\([,:[:space:]]\\)|\\1FL_RD\\2|' \\
                -e 's|\\([=:[:space:]]\\)AD\\([,:[:space:]]\\)|\\1FL_AD\\2|' \\
                -e 's|\\([=:[:space:]]\\)VF\\([,:[:space:]]\\)|\\1FL_VF\\2|' \\
                -e 's|\\([=:[:space:]]\\)DPP\\([,:[:space:]]\\)|\\1FL_DPP\\2|' \\
                -e 's|\\([=:[:space:]]\\)DPN\\([,:[:space:]]\\)|\\1FL_DPN\\2|' \\
                -e 's|\\([=:[:space:]]\\)RDP\\([,:[:space:]]\\)|\\1FL_RDP\\2|' \\
                -e 's|\\([=:[:space:]]\\)RDN\\([,:[:space:]]\\)|\\1FL_RDN\\2|' \\
                -e 's|\\([=:[:space:]]\\)ADP\\([,:[:space:]]\\)|\\1FL_ADP\\2|' \\
                -e 's|\\([=:[:space:]]\\)ADN\\([,:[:space:]]\\)|\\1FL_ADN\\2|' \\
                "\${fillout_vcf}" > "\${fillout_relabel_vcf}"
                bgzip -c "\${fillout_relabel_vcf}" > "\${fillout_relabel_vcf_gz}"
                tabix "\${fillout_relabel_vcf_gz}"
                # combine the two sets of vcfs
                # pull off the merge vcf header and create an annotations file
                merged_vcf="${ return inputs.merged_vcf.path }"
                merged_vcf_header="merged.header.txt"
                merged_vcf_gz="${ return inputs.merged_vcf_gz.path }"
                fillout_merged_vcf="fillout.merged.vcf"
                grep -E '##FORMAT|##INFO' "\${merged_vcf}" > "\${merged_vcf_header}"
                bcftools annotate --header-lines "\${merged_vcf_header}" --annotations "\${merged_vcf_gz}" --columns 'FORMAT/GT,FORMAT/AD,FORMAT/DP,INFO/AC,INFO/AN' --output "\${fillout_merged_vcf}" --output-type v "\${fillout_relabel_vcf_gz}"
                # label each variant with its source samples; samples with '.' empty AD value weren't present in sample originally
                # create a new annotations file that has the samples listed for each variant in a new INFO field labeled SRC
                fillout_merged_vcf_header="fillout.merged.header.txt"
                fillout_merged_annotation="fillout.merged.annotation.txt"
                fillout_merged_annotation_gz="fillout.merged.annotation.txt.gz"

                #  ! ! ! ! ! IMPORTANT ! ! ! ! !
                fillout_merged_samplesources="fillout.merged.sources.vcf"
                #  ^^^ This is going to be the final output file at this step
                #  containing the merge product of the fillout vcf + original sample vcf's + SRC sample source labels per-variant
                #  ! ! ! ! ! ! ! ! ! !

                # get the old header lines
                grep '##' "\${fillout_merged_vcf}" > "\${fillout_merged_vcf_header}"
                # add a new line
                echo '##INFO=<ID=SRC,Type=String,Number=.,Description="Source samples for the variant">' >> "\${fillout_merged_vcf_header}"
                # start the annotation file header
                echo '#CHROM\tPOS\tREF\tAL\tINFO' > "\${fillout_merged_annotation}"
                # get all the variant sample labels for variants with AD>0;
                # this means they were from the sample originally and not from fillout
                bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%SAMPLE,]\\n' -i 'AD>0' "\${fillout_merged_vcf}" >> "\${fillout_merged_annotation}"
                bgzip -c "\${fillout_merged_annotation}" > "\${fillout_merged_annotation_gz}"
                tabix -p vcf "\${fillout_merged_annotation_gz}"
                # re-annotate the final file with the new header lines from the fillout file
                bcftools annotate --header-lines "\${fillout_merged_vcf_header}" \\
                --annotations "\${fillout_merged_annotation_gz}" \\
                -c CHROM,POS,REF,ALT,SRC \\
                --output "\${fillout_merged_samplesources}" \\
                --output-type v \\
                "\${fillout_merged_vcf}"
      inputs:
        fillout_vcf: File # output from GetBaseCountsMultiSample
        merged_vcf: File # the vcf merged from all sample inputs # TODO: do not actually need this could just use the gz version
        merged_vcf_gz:
          type: File
          secondaryFiles:
            - .tbi
      outputs:
        fillout_sources_vcf:
          type: File
          outputBinding:
            glob: fillout.merged.sources.vcf


  # FILLOUT POST-PROCESSING STARTS HERE
  fillout_post_processing:
    run: fillout_post_processing.cwl
    in:
      samples: samples
      fillout_vcf: fix_labels_and_merge_vcfs/fillout_sources_vcf
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
    outputSource: fillout_post_processing/portal_file # called_file
  uncalled_file:
    type: File
    outputSource: fillout_post_processing/uncalled_file
