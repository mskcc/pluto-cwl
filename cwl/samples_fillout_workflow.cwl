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
  samples: "types.yml#FilloutMafOptionalSample[]"
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
  # Check for samples that lack a .maf file
  get_samples_with_without_maf:
    doc: separate the samples list into samples with and without .maf files
    in:
      samples: samples
    out: [ samples_with_maf, samples_without_maf ]
    run: fillout_get_samples_with_without_maf.cwl

  get_bam_and_ids_from_samples_without_maf:
    doc:
    in:
      samples: get_samples_with_without_maf/samples_without_maf
    out: [ sample_ids, bam_files ]
    run:
      class: ExpressionTool
      inputs:
        samples:
          type: "types.yml#FilloutNoMafsample[]"
      outputs:
        sample_ids: string[]
        bam_files: File[]
      expression: |
        ${
          var sample_ids = [];
          var bam_files = [];

          for ( var i in inputs.samples ){
            sample_ids.push(inputs.samples[i]['sample_id']);
            bam_files.push(inputs.samples[i]['bam_file']);
          };

          var res = {
              "sample_ids": sample_ids,
              "bam_files": bam_files
            };

          return res;
        }


  # PRE-PROCESSING
  fillout_pre_processing:
    run: fillout_pre_processing.cwl
    in:
      samples: get_samples_with_without_maf/samples_with_maf
      ref_fasta: ref_fasta
    out:
      [ sample_ids, clinical_sample_ids, bam_files, vcf_gz_files, merged_vcf, merged_vcf_gz ]


  # PRIMARY FILLOUT PROCESSING STARTS HERE
  # TODO: convert this to a `scatter` step that runs per-sample in parallel, then merge the outputs otherwise we will hit the command line arg length issues eventually
  # NOTE: maybe do not do this ^^^ because currently its useful to have a multi-sample vcf output, need to investigate this more
  gbcms:
    doc: run GetBaseCountsMultiSample on all the bam files against the target regions (the merged vcf from all samples)
    in:
      sample_ids:
        source: [ fillout_pre_processing/sample_ids,  get_bam_and_ids_from_samples_without_maf/sample_ids ]
        linkMerge: merge_flattened
      bam_files:
        source: [ fillout_pre_processing/bam_files, get_bam_and_ids_from_samples_without_maf/bam_files ]
        linkMerge: merge_flattened
      targets_vcf: fillout_pre_processing/merged_vcf
      ref_fasta: ref_fasta
    out: [ output_file ]
    run:
      class: CommandLineTool
      id: gbcms
      label: gbcms
      baseCommand: [ "sh", "run.gbcms.sh" ]
      requirements:
        InlineJavascriptRequirement: {}
        ResourceRequirement:
          coresMin: 8
        DockerRequirement:
          dockerPull: cmopipeline/getbasecountsmultisample:1.2.2
        InitialWorkDirRequirement:
          listing:
            - entryname: run.gbcms.sh
              entry: |-
                set -eu
                # make the bam arg string; looks like this:
                # --bam Sample5:Sample5.rg.md.abra.printreads.bam --bam Sample6:Sample6.rg.md.abra.printreads.bam
                bams_arg='${
                  var args1 = inputs.sample_ids.map((k, i) => [`${k}:${inputs.bam_files[i].path}`])
                  var args2 = args1.map( (a) => "--bam " + a )
                  return args2.join(" ") ;
                  }'
                sample_ids="${ return inputs.sample_ids.join("\t"); }"
                fasta="${ return inputs.ref_fasta.path; }"
                vcf="${ return inputs.targets_vcf.path }"
                fillout_vcf="fillout.vcf"
                num_lines="\$(grep -v '^[#]' \${input_maf} | wc -l)"
                more_than_zero="\$(( \${num_lines} > 0))"
                if [ "\${more_than_zero}" == "1" ]
                then
                  GetBaseCountsMultiSample --fasta "\${fasta}" --vcf "\${vcf}" --maq 20 --baq 20 --filter_improper_pair 0 --thread 8 --output "\${fillout_vcf}" \${bams_arg}
                else
                  cat << EOF > "\${fillout_vcf}"
                ##fileformat=VCFv4.2
                ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
                ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth matching reference (REF) allele">
                ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth matching alternate (ALT) allele">
                ##FORMAT=<ID=VF,Number=1,Type=Float,Description="Variant frequence (AD/DP)">
                ##FORMAT=<ID=DPP,Number=1,Type=Integer,Description="Depth on postitive strand">
                ##FORMAT=<ID=DPN,Number=1,Type=Integer,Description="Depth on negative strand">
                ##FORMAT=<ID=RDP,Number=1,Type=Integer,Description="Reference depth on postitive strand">
                ##FORMAT=<ID=RDN,Number=1,Type=Integer,Description="Reference depth on negative strand">
                ##FORMAT=<ID=ADP,Number=1,Type=Integer,Description="Alternate depth on postitive strand">
                ##FORMAT=<ID=ADN,Number=1,Type=Integer,Description="Alternate depth on negative strand">
                ##FORMAT=<ID=DPF,Number=1,Type=Integer,Description="Total fragment depth">
                ##FORMAT=<ID=RDF,Number=1,Type=Float,Description="Fragment depth matching reference (REF) allele">
                ##FORMAT=<ID=ADF,Number=1,Type=Float,Description="Fragment depth matching alternate (ALT) allele">
                #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  P-\${sample_ids}
                EOF
                fi
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
    # NOTE: IMPORTANT: Need to apply FL_VF, and SRC fields !! Those are the most imporant fields for downstream uses
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
      baseCommand: [ "bash", "run.fix_labels_and_merge_vcfs.sh" ]
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:21.4.1
        InitialWorkDirRequirement:
          listing:
            - entryname: run.fix_labels_and_merge_vcfs.sh
              entry: |-
                set -eu
                # need to rename some conflicting sample INFO tags; prepend FL for "fillout"
                # change these tags: DP:RD:AD:VF:DPP:DPN:RDP:RDN:ADP:ADN
                # pre-pend them all with 'FL_' for fillout so they do not collide with downstream labels

                # input file: fillout.vcf
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

                # input file: merged.vcf
                merged_vcf="${ return inputs.merged_vcf.path }"
                # should contain these lines;
                ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths of REF and ALT(s) in the order listed">
                ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
                ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
                ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">

                # singleton Sample2.sorted.vcf has these lines;
                ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths of REF and ALT(s) in the order listed">
                ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">

                merged_vcf_header="merged.header.txt"
                merged_vcf_gz="${ return inputs.merged_vcf_gz.path }"
                fillout_merged_vcf="fillout.merged.vcf"

                # pull off the input vcf header for annotating
                grep -E '##FORMAT|##INFO' "\${merged_vcf}" > "\${merged_vcf_header}"

                # need to check that these header lines exist because they tend to not exist if the input was a singleton sample
                if ! grep -q ID=AC "\${merged_vcf_header}"; then
                echo '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">' >> "\${merged_vcf_header}"
                fi

                if ! grep -q ID=AN "\${merged_vcf_header}"; then
                echo '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">' >> "\${merged_vcf_header}"
                fi

                # apply the annotations from the merged vcf to the fillout vcf
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
                # NOTE: 'ALT' column is misspelled here!! Need to fix this!!
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
        fillout_sources_vcf: # this is used as the fillout_vcf downstream
          type: File
          outputBinding:
            glob: fillout.merged.sources.vcf



  # # FILLOUT POST-PROCESSING STARTS HERE
  # we do not need the input .maf files anymore so convert all samples to FilloutNoMafsample format
  convert_all_to_FilloutNoMafsample:
    in:
      samples_with_maf: get_samples_with_without_maf/samples_with_maf
      samples_without_maf: get_samples_with_without_maf/samples_without_maf
    out: [ samples ]
    run:
      class: ExpressionTool
      inputs:
        samples_with_maf: "types.yml#FilloutSample[]"
        samples_without_maf: "types.yml#FilloutNoMafsample[]"
      outputs:
        samples: "types.yml#FilloutNoMafsample[]"
      expression: |
        ${
          var samples = [];

          for ( var i in inputs.samples_without_maf ){
            samples.push(inputs.samples_without_maf[i]);
          };

          for ( var i in inputs.samples_with_maf ){
            var sample = {
              "sample_id": inputs.samples_with_maf[i]["sample_id"],
              "normal_id": inputs.samples_with_maf[i]["normal_id"],
              "sample_type": inputs.samples_with_maf[i]["sample_type"],
              "bam_file": inputs.samples_with_maf[i]["bam_file"],
            };
            samples.push(sample);
          };

          var res = {
              "samples": samples,
            };

          return res;
        }

  fillout_post_processing:
    doc: run post processing on the files based on the metadata bundled in with the samples
    run: fillout_post_processing.cwl
    in:
      samples: convert_all_to_FilloutNoMafsample/samples
      fillout_vcf: fix_labels_and_merge_vcfs/fillout_sources_vcf
      clinical_sample_ids: fillout_pre_processing/clinical_sample_ids
      ref_fasta: ref_fasta
      exac_filter: exac_filter
    out: [ output_file, filtered_file, portal_file, uncalled_file ]





outputs:
  merged_vcf:
    doc: all the input maf files merged together and converted back to vcf with duplicate variants removed
    type: File
    outputSource: fillout_pre_processing/merged_vcf
  fillout_sources_vcf:
    doc: the GetBaseCountsMultiSample output file with the SRC originating sample labels added for each variant entry
    type: File
    outputSource: fix_labels_and_merge_vcfs/fillout_sources_vcf
  output_file:
    doc: the primary output maf file with all filled out variants
    type: File
    outputSource: fillout_post_processing/output_file
  filtered_file:
    doc: the output file after applying some filter criteria
    type: File
    outputSource: fillout_post_processing/filtered_file
  portal_file:
    doc: the filtered file after applying cBioPortal formatting
    type: File
    outputSource: fillout_post_processing/portal_file # called_file
  uncalled_file:
    doc: the mutations that were removed from the portal file
    type: File
    outputSource: fillout_post_processing/uncalled_file
