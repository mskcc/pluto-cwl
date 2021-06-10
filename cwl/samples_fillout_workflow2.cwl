#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
doc: "
Workflow to run GetBaseCountsMultiSample fillout on a number of samples, each with their own bam and maf files
"
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  # NOTE: arrays for sample_ids, bam_files, maf_files must all be the same length and in the same order by sample
  sample_ids:
    type:
        type: array
        items: string
  bam_files:
    type:
        type: array
        items: File
    secondaryFiles:
        - ^.bai
  maf_files:
    type:
        type: array
        items: File
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

steps:
  # convert all maf input files back to vcf because they are much easier to manipulate that way
  # NOTE: This is important; do NOT try to do these complex manipulations on maf format file
  maf2vcf:
    scatter: [ sample_id, maf_file ]
    scatterMethod: dotproduct
    in:
      sample_id: sample_ids
      maf_file: maf_files
      ref_fasta: ref_fasta
    out:
      [ output_file ]
    run:
      class: CommandLineTool
      baseCommand: ['bash', 'run.sh']
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:latest
        InitialWorkDirRequirement:
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

  # merge all the vcf files together to create a vcf with all of the variant regions for all samples
  # this will be used as the target regions for fillout
  merge_vcfs:
    in:
      sample_ids: sample_ids
      vcf_gz_files: maf2vcf/output_file
    out:
      [ merged_vcf, merged_vcf_gz ]
    run:
      class: CommandLineTool
      baseCommand: ['bash', 'run.sh']
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:latest
        InitialWorkDirRequirement:
          listing:
            - entryname: run.sh
              entry: |-
                set -eu
                # get a comma-delim string of the sample names
                samples_arg="${ return inputs.sample_ids.join(',') }"
                # get a space-delim string of vcf file paths
                vcf_files="${ return inputs.vcf_gz_files.map((a) => a.path).join(' ') }"
                merged_vcf="merged.vcf"
                merged_vcf_gz="merged.vcf.gz"
                bcftools merge --force-samples --output-type v \${vcf_files} | \\
                bcftools sort | \\
                bcftools view -s "\${samples_arg}" > "\${merged_vcf}"
                bgzip -c "\${merged_vcf}" > "\${merged_vcf_gz}"
                tabix "\${merged_vcf_gz}"
      inputs:
        sample_ids: string[]
        vcf_gz_files:
          type:
              type: array
              items: File
          secondaryFiles:
              - .tbi
      outputs:
        merged_vcf:
          type: File
          outputBinding:
            glob: merged.vcf
        merged_vcf_gz:
          type: File
          outputBinding:
            glob: merged.vcf.gz
          secondaryFiles:
            - .tbi

  # run GetBaseCountsMultiSample on all the bam files against the target regions (the merged vcf from all samples)
  gbcms:
    in:
      sample_ids: sample_ids
      bam_files: bam_files
      targets_vcf: merge_vcfs/merged_vcf
      ref_fasta: ref_fasta
    out: [ output_file ]
    run:
      class: CommandLineTool
      baseCommand: [ "sh", "run.sh" ]
      requirements:
        InlineJavascriptRequirement: {}
        ResourceRequirement:
          coresMin: 4
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
                GetBaseCountsMultiSample --fasta "\${fasta}" --vcf "\${vcf}" --maq 20 --baq 20 --filter_improper_pair 0 --thread 4 --output "\${fillout_vcf}" \${bams_arg}
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

  # need to perform a bunch of steps to clean up the fillout vcf and re-annotate it with sample labels
  # also we are going to merge in the "merged_vcf" with the original samples' vcf values
  # into the fillout vcf so we have both old and new (sample + fillout) values in a single vcf
  # also we are going to add a column called SRC telling the source (which sample) each variant was originally found in
  fix_labels_and_merge_vcfs:
    in:
      sample_ids: sample_ids
      fillout_vcf: gbcms/output_file
      merged_vcf: merge_vcfs/merged_vcf
      merged_vcf_gz: merge_vcfs/merged_vcf_gz
      ref_fasta: ref_fasta
    out: [ fillout_sources_vcf ]
    run:
      class: CommandLineTool
      baseCommand: [ "bash", "run.sh" ]
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:latest
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
                --output-type v "\${fillout_merged_vcf}"
              # # convert the multi-sample annotated fillout vcf back into individual sample maf files
              # for sampleid in \${samples_arg}; do
              # sample_maf="\${sampleid}.fillout.maf"
              # vcf2maf.pl \\
              # --input-vcf "\${fillout_merged_samplesources}" \\
              # --output-maf "\${sample_maf}" \\
              # --ref-fasta "\${ref_fasta}" \\
              # --cache-version 86 \\
              # --species homo_sapiens \\
              # --ncbi-build GRCh37 \\
              # --vep-path /usr/bin/vep \\
              # --vep-data /var/cache \\
              # --filter-vcf /juno/work/ci/resources/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \\
              # --retain-info AC,AN,SRC \\
              # --retain-fmt GT,FL_AD,FL_ADN,FL_ADP,FL_DP,FL_DPN,FL_DPP,FL_RD,FL_RDN,FL_RDP,FL_VF,AD,DP \\
              # --vcf-tumor-id \\
              # "\${sampleid}" \\
              # --tumor-id "\${sampleid}"
              # done
      inputs:
        fillout_vcf: File # output from GetBaseCountsMultiSample
        merged_vcf: File # the vcf merged from all sample inputs # TODO: do not actually need this could just use the gz version
        merged_vcf_gz:
          type: File
          secondaryFiles:
            - .tbi
        # sample_ids: string[]
        # ref_fasta:
        #   type: File
        #   secondaryFiles:
        #     - .amb
        #     - .ann
        #     - .bwt
        #     - .pac
        #     - .sa
        #     - .fai
        #     - ^.dict
        # get a space-delim string of the sample names
        # samples_arg="${ return inputs.sample_ids.join(' ') }"
        # ref_fasta="${ return inputs.ref_fasta.path ; }"
      outputs:
        fillout_sources_vcf:
          type: File
          outputBinding:
            glob: fillout.merged.sources.vcf

  # next we need to split apart the merged fillout vcf back into individual sample maf files
  # split_merged_vcf:
  #   scatter: [ sample_id ]
  #   in:
  #     sample_id: sample_ids
  #     fillout_vcf: fix_labels_and_merge_vcfs/fillout_sources_vcf
  #   out: []
  #   run: []





outputs:
  fillout_vcf:
    type: File
    outputSource: fix_labels_and_merge_vcfs/fillout_sources_vcf
