#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
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
  # create a list of just sample_ids out of the samples record array
  create_samples_list:
    in:
      samples: samples
    out: [ sample_ids ]
    run:
      class: ExpressionTool
      inputs:
        samples: "types.yml#FilloutSample[]"
      outputs:
        sample_ids: string[]
      # NOTE: in the line below `var i in inputs.samples`, `i` is an int representing the index position in the array `inputs.samples`
      # in Python it would look like ` x = ['a', 'b']; for i in range(len(x)): print(i, x[i]) `
      expression: "${
        var sample_ids = [];
        for ( var i in inputs.samples ){
            sample_ids.push(inputs.samples[i]['sample_id']);
          };
        return {'sample_ids': sample_ids};
        }"

  # get a list of just the bam files for downstream process
  create_bam_list:
    in:
      samples: samples
    out:
      [ bam_files ]
    run:
      class: ExpressionTool
      inputs:
        samples: "types.yml#FilloutSample[]"
      outputs:
        bam_files:
          type:
              type: array
              items: File
          secondaryFiles:
              - ^.bai
      expression: "${
        var bam_files = [];
        for ( var i in inputs.samples ){
            bam_files.push(inputs.samples[i]['bam_file']);
          };
        return {'bam_files': bam_files};
        }"


  # convert all maf input files back to vcf because they are much easier to manipulate that way
  # NOTE: This is important; do NOT try to do complex manipulations on maf format file
  maf2vcf:
    scatter: sample
    in:
      sample: samples
      sample_id:
        valueFrom: ${ return inputs.sample['sample_id']; }
      maf_file:
        valueFrom: ${ return inputs.sample['maf_file']; }
      ref_fasta: ref_fasta
    out:
      [ output_file ]
    run:
      class: CommandLineTool
      baseCommand: ['bash', 'run.sh']
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:21.4.1
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
  # this will be used as the target regions for fillout (GetBaseCountsMultiSample)
  merge_vcfs:
    in:
      sample_ids: create_samples_list/sample_ids
      vcf_gz_files: maf2vcf/output_file
    out:
      [ merged_vcf, merged_vcf_gz ]
    run:
      class: CommandLineTool
      baseCommand: ['bash', 'run.sh']
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:21.4.1
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
  # TODO: convert this to a `scatter` step that runs per-sample in parallel, then merge the outputs
  # otherwise we will hit the command line arg length issues
  gbcms:
    in:
      sample_ids: create_samples_list/sample_ids
      bam_files: create_bam_list/bam_files
      targets_vcf: merge_vcfs/merged_vcf
      ref_fasta: ref_fasta
    out: [ output_file ]
    run:
      class: CommandLineTool
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

  # need to perform a bunch of steps to clean up the fillout vcf and re-annotate it with sample labels
  # also we are going to merge in the "merged_vcf" with the original samples' vcf values
  # into the fillout vcf so we have both old and new (sample + fillout) values in a single vcf
  # also we are going to add a column called SRC telling the source (which sample) each variant was originally found in
  fix_labels_and_merge_vcfs:
    in:
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


  # create CLI args for bcftools filter in order to filter on sample ID's in SRC INFO tag
  # use the bash snippet generated here in a downstream process for filtering
  # should look like this:
  # 'SRC="Sample1,Sample2,"'
  make_filter_args:
    in:
      samples: samples
    out: [ research_arg, clinical_arg, clinical_expression ]
    run:
      class: ExpressionTool
      inputs:
        samples: "types.yml#FilloutSample[]"
      outputs:
        # 'SRC="Sample1,Sample2,"'
        research_arg: string
        # 'SRC="Sample1,Sample2,"'
        clinical_arg: string
        # bcftools view -e "SRC='Sample3'" - |  bcftools view -e "SRC='Sample4'" - | ...
        clinical_expression: string
      expression: |
        ${
          // split the list of samples into research and clinical samples
          var research_samples = [];
          var clinical_samples = [];

          for ( var i in inputs.samples ){

            var sample_type = inputs.samples[i]['sample_type'];
            var sample_id = inputs.samples[i]['sample_id'];

            if ( sample_type == 'research' ){
              research_samples.push(sample_id);

            } else if ( sample_type == 'clinical' ) {
              clinical_samples.push(sample_id);
            }

          };

          // create the bcftools filter expression for samples;
          // SRC='Sample1,Sample2'
          var research_arg = `SRC='${research_samples.join(',')}'`;
          var clinical_arg = `SRC='${clinical_samples.join(',')}'`;





          // create a set of bcftools view -e commands for the clinical samples
          // need to chain together individual exclusions for each clinical sample_id
          // it will look like
          // bcftools view -e "SRC='Sample3'" - |  bcftools view -e "SRC='Sample4'" - | ...

          // start with empty string and build up the final bash command
          var expr = '';

          // check if there are >0
          var num_clinical_samples = clinical_samples.length;
          if (num_clinical_samples > 0){
            // need to start the command with a | character if there are any samples
            expr = expr + ' | '
          };

          for ( var i in clinical_samples ){
            var sample_id = clinical_samples[i];

            // be careful with the quoting here; must be "SRC='Sample1'"
            expr = expr + 'bcftools view -e "SRC=' + "'" + sample_id +"'" + '" ';

            // need to apply the | between all args except the final one
            if (i < num_clinical_samples - 1){
              expr = expr + ' | ';
            };
          };

          return {'research_arg': research_arg, 'clinical_expression': expr, 'clinical_arg': clinical_arg};
        }

  # split the multi-sample vcf into individual per-sample vcf files and apply conditional filters
  # TODO: deprecate the need for bcftools +split by making GetBaseCountsMultiSample run individually per-sample
  split_filter_vcf:
    scatter: sample
    in:
      sample: samples
      fillout_vcf: fix_labels_and_merge_vcfs/fillout_sources_vcf
      research_arg: make_filter_args/research_arg
      clinical_arg: make_filter_args/clinical_arg
      clinical_expression: make_filter_args/clinical_expression
    out: [ sample ]
    run:
      class: CommandLineTool
      baseCommand: [ "bash", "run.sh" ]
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:21.4.1
        InitialWorkDirRequirement:
          listing:
            - $(inputs.fillout_vcf)
            - entryname: run.sh
              entry: |-
                set -eux

                input_vcf='${ return inputs.fillout_vcf.basename ; }'
                sample_type='${ return inputs.sample["sample_type"] ; }'
                sample_id='${ return inputs.sample["sample_id"] ; }'
                # be careful with the quoting args here;
                # research_filter_arg="${ return inputs.research_arg ; }"
                # clinical_expression='${return inputs.clinical_expression ; }'

                # dir to hold split vcf contents
                split_dir="split"
                split_vcf="\${split_dir}/\${sample_id}.vcf"

                # file where we will save the final filtered variants for conversion to maf
                filtered_vcf="\${sample_id}.filtered.vcf"
                # file where we will write out variants to exclude based on filter criteria
                filter_exclusion_file="\${sample_id}.exclude.txt"
                # file to save a .txt version of the .vcf for intersecting with the exclusion list
                sample_vcf_txt="\${sample_id}.vcf.txt"
                # file where we will save the list of variants to keep
                filter_inclusion_file="\${sample_id}.include.txt"

                # final output filename
                unfiltered_vcf="\${sample_id}.vcf"

                mkdir -p "\${split_dir}"

                # split the multi-sample vcf into per-sample vcf files
                bcftools +split "\${input_vcf}" --output-type v --output "\${split_dir}"

                # CONDITIONAL VARIANT FILTERING
                # in clinical samples ONLY;
                # exclude filled-out variants (identified by invalid AD value of ".")
                # from research samples (not in any clinical samples)
                # that have VAF >0.1
                if [ "\${sample_type}" == "clinical" ]; then


                # NOTE: Filter feature is not ready yet; just output the same file twice for now
                # # keep this section as a placeholder for the upcoming filter feature
                # # keep these code blocks as examples of how to implement
                cp "\${split_vcf}" "\${unfiltered_vcf}"
                cp "\${unfiltered_vcf}" "\${filtered_vcf}"

                # # NEW METHOD
                # bcftools filter -e "${ return inputs.clinical_arg; } && FL_VF<0.1 & AD='.'" "\${split_vcf}" > "\${filtered_vcf}"
                # # NOTE: add a bcftools query here for an extra .tsv output file for convenience
                #
                #
                #
                # # OLD METHOD:
                # #
                # # # create file to hold list of variants to exclude
                # # # exclusion list criteria:
                # # # variants with fillout VAF >0.1
                # # # that were not present in the original sample (invalid 'AD' value = "is_fillout=TRUE")
                # # # which originated in some research sample
                # # # but were not present in any clinical sample ('SRC' sample list)
                # # # # NOTE: https://samtools.github.io/bcftools/bcftools.html#expressions
                # # # # > Comma in strings is interpreted as a separator and when multiple values are compared, the OR logic is used
                # # bcftools view -i 'FL_VF>0.1' "\${split_vcf}" | \
                # # bcftools view -i 'AD="."' | \
                # # bcftools view -i "${ return inputs.research_arg ; }" ${return inputs.clinical_expression ; } | \
                # # bcftools query -f '%CHROM\\t%POS\\t%END\\t%SRC\\t[AD=%AD\\tFL_VF=%FL_VF]\\n' - > "\${filter_exclusion_file}"
                # #
                # # # convert the original vcf to a flat txt format for use with bedtools intersect
                # # bcftools query -f '%CHROM\\t%POS\\t%END\\t%SRC\\t[AD=%AD\\tFL_VF=%FL_VF]\\n' "\${split_vcf}" > "\${sample_vcf_txt}"
                # #
                # # # get the list of variants that are not in the exclusion list
                # # bedtools intersect -v -a "\${sample_vcf_txt}" -b "\${filter_exclusion_file}" -wa > "\${filter_inclusion_file}"
                # #
                # # # filter the original vcf down to only include the loci from the inclusion file
                # # bcftools filter --targets-file "\${filter_inclusion_file}" "\${split_vcf}" > "\${filtered_vcf}"
                #
                # # move the unfiltered_vcf to the designated output path
                # cp "\${split_vcf}" "\${unfiltered_vcf}"

                else

                # sample did not need filtering; copy the sample's vcf to the designated path
                cp "\${split_vcf}" "\${unfiltered_vcf}"
                cp "\${unfiltered_vcf}" "\${filtered_vcf}"

                fi
      inputs:
        sample: "types.yml#FilloutSample"
        clinical_expression: string
        research_arg: string
        clinical_arg: string
        # multi-sample vcf generated from GetBaseCountsMultiSample on all samples at once
        fillout_vcf: File
      outputs:
        unfiltered_vcf:
          type: File
          outputBinding:
            glob: $(inputs.sample['sample_id']).vcf
        filtered_vcf:
          type: File
          outputBinding:
            glob: $(inputs.sample['sample_id']).filtered.vcf
        sample:
          type: "types.yml#FilloutSample"
          outputBinding:
            # set the unfiltered_vcf to the unfiltered_vcf
            outputEval: ${
              var ret = inputs.sample;
              ret['unfiltered_vcf'] = {"class":"File", "path":runtime.outdir + "/" + inputs.sample["sample_id"] + ".vcf"};
              ret['filtered_vcf'] = {"class":"File", "path":runtime.outdir + "/" + inputs.sample["sample_id"] + ".filtered.vcf"};
              return ret;
              }


  # convert each sample filtered vcf back to maf format for cBioPortal and end user
  vcf_to_maf:
    scatter: sample
    in:
      sample: split_filter_vcf/sample
      ref_fasta: ref_fasta
      exac_filter: exac_filter
    out: [ unfiltered_maf, filtered_maf, sample ]
    run:
      class: CommandLineTool
      baseCommand: [ "sh", "run.sh" ]
      requirements:
        ResourceRequirement:
          coresMin: 8
        DockerRequirement:
          dockerPull: mskcc/roslin-variant-vcf2maf:1.6.17
        InitialWorkDirRequirement:
          listing:
            - $(inputs.sample['unfiltered_vcf'])
            - $(inputs.sample['filtered_vcf'])
            - entryname: run.sh
              entry: |-
                set -eu
                # convert the multi-sample annotated fillout vcf back into individual sample maf files
                sample_id="${ return inputs.sample['sample_id'] ; }"
                normal_id="${ return inputs.sample['normal_id'] ; }"
                ref_fasta="${ return inputs.ref_fasta.path ; }"
                unfiltered_vcf="${ return inputs.sample['unfiltered_vcf'].basename ; }"
                filtered_vcf="${ return inputs.sample['filtered_vcf'].basename ; }"
                exac_filter="${ return inputs.exac_filter.path ; }"

                # not sure why I have to do this but if I dont then it breaks looking for some file;
                # ERROR: Couldn't open VCF: /var/lib/cwl/stg640f09da-4a21-46eb-bc90-a27199b424bc/fillout.merged.sources.split.vcf!
                # cp "\${unfiltered_vcf}" input.vcf

                # main output files:
                unfiltered_maf="\${sample_id}.fillout.maf"
                filtered_maf="\${sample_id}.fillout.filtered.maf"

                # NOTE: /usr/bin/vep, /var/cache is inside the container already
                perl /usr/bin/vcf2maf/vcf2maf.pl \\
                --input-vcf "\${unfiltered_vcf}" \\
                --output-maf "\${unfiltered_maf}" \\
                --ref-fasta "\${ref_fasta}" \\
                --cache-version 86 \\
                --species homo_sapiens \\
                --ncbi-build GRCh37 \\
                --vep-path /usr/bin/vep \\
                --vep-data /var/cache \\
                --filter-vcf "\${exac_filter}" \\
                --retain-info AC,AN,SRC \\
                --retain-fmt GT,FL_AD,FL_ADN,FL_ADP,FL_DP,FL_DPN,FL_DPP,FL_RD,FL_RDN,FL_RDP,FL_VF,AD,DP \\
                --vep-forks 8 \\
                --vcf-tumor-id "\${sample_id}" \\
                --tumor-id "\${sample_id}" \\
                --normal-id "\${normal_id}"

                perl /usr/bin/vcf2maf/vcf2maf.pl \\
                --input-vcf "\${filtered_vcf}" \\
                --output-maf "\${filtered_maf}" \\
                --ref-fasta "\${ref_fasta}" \\
                --cache-version 86 \\
                --species homo_sapiens \\
                --ncbi-build GRCh37 \\
                --vep-path /usr/bin/vep \\
                --vep-data /var/cache \\
                --filter-vcf "\${exac_filter}" \\
                --retain-info AC,AN,SRC \\
                --retain-fmt GT,FL_AD,FL_ADN,FL_ADP,FL_DP,FL_DPN,FL_DPP,FL_RD,FL_RDN,FL_RDP,FL_VF,AD,DP \\
                --vep-forks 8 \\
                --vcf-tumor-id "\${sample_id}" \\
                --tumor-id "\${sample_id}" \\
                --normal-id "\${normal_id}"
      inputs:
        sample: "types.yml#FilloutSample"
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
        exac_filter: File
        # NOTE: this might need secondaryFiles
        # Could not load .tbi/.csi index of /var/lib/cwl/stg2e0eeed1-680c-4c33-beeb-c4dd90309998/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
      outputs:
        unfiltered_maf:
          type: File
          outputBinding:
            glob: $(inputs.sample['sample_id']).fillout.maf
        filtered_maf:
          type: File
          outputBinding:
            glob: $(inputs.sample['sample_id']).fillout.filtered.maf
        sample:
          type: "types.yml#FilloutSample"
          outputBinding:
            # set the unfiltered_maf for the sample
            outputEval: ${
              var ret = inputs.sample;
              ret['unfiltered_maf'] = {"class":"File", "path":runtime.outdir + "/" + inputs.sample["sample_id"] + ".fillout.maf"};
              ret['filtered_maf'] = {"class":"File", "path":runtime.outdir + "/" + inputs.sample["sample_id"] + ".fillout.filtered.maf"};
              return ret;
              }

  # combine all the individual mafs into a single maf; add comment headers; fix some values
  concat_with_comments:
    in:
      unfiltered_mafs: vcf_to_maf/unfiltered_maf
      filtered_mafs: vcf_to_maf/filtered_maf
      output_unfiltered_filename: output_fname
    out: [ output_file, filtered_file ]
    run:
      class: CommandLineTool
      baseCommand: [ "bash", "run.sh" ]
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:21.4.1
        InitialWorkDirRequirement:
          listing:
            - entryname: run.sh
              entry: |-
                set -eu
                # get a space-delim string of file paths
                input_unfiltered_files="${ return inputs.unfiltered_mafs.map((a) => a.path).join(' ') }"
                input_filtered_files="${ return inputs.filtered_mafs.map((a) => a.path).join(' ') }"
                output_unfiltered_filename="${ return inputs.output_unfiltered_filename }"
                output_filtered_filename="${ return inputs.output_filtered_filename }"

                concat_fillout_maf=fillout.maf
                concat_filtered_maf=fillout.filtered.maf

                concat-tables.py -o "\${concat_fillout_maf}" --no-carriage-returns --comments --progress --na-str '' \${input_unfiltered_files}

                concat-tables.py -o "\${concat_filtered_maf}" --no-carriage-returns --comments --progress --na-str '' \${input_filtered_files}

                # fix issues with blank ref alt count cols for fillout variants
                fillout_tmp=tmp.tsv
                fillout_filtered_tmp=tmp.filtered.tsv

                update_fillout_maf.py "\${concat_fillout_maf}" "\${fillout_tmp}"
                update_fillout_maf.py "\${concat_filtered_maf}" "\${fillout_filtered_tmp}"

                # need to add header comments for new columns
                echo '# SRC="Samples variant was originally found in (Fillout)"' > comments
                echo '# FL_AD="Depth matching alternate (ALT) allele (Fillout)"' >> comments
                echo '# FL_ADN="Alternate depth on negative strand (Fillout)"' >> comments
                echo '# FL_ADP="Alternate depth on postitive strand (Fillout)"' >> comments
                echo '# FL_DP="Total depth (Fillout)"' >> comments
                echo '# FL_DPN="Depth on negative strand (Fillout)"' >> comments
                echo '# FL_DPP="Depth on postitive strand (Fillout)"' >> comments
                echo '# FL_RD="Depth matching reference (REF) allele (Fillout)"' >> comments
                echo '# FL_RDN="Reference depth on negative strand (Fillout)"' >> comments
                echo '# FL_RDP="Reference depth on postitive strand (Fillout)"' >> comments
                echo '# FL_VF="Variant frequence (AD/DP) (Fillout)"' >> comments
                echo '# AD="Allelic Depths of REF and ALT(s) in the order listed (Sample)"' >> comments
                echo '# DP="Read Depth (Sample)"' >> comments
                echo '# depth="Final Read Depth value; Sample Depth if present, otherwise based on Fillout Depth"' >> comments
                echo '# ref_count="Final Allelic Depth of REF; Sample Allelic Depth if present, otherwise based on Fillout Allelic Depth"' >> comments
                echo '# alt_count="Final Allelic Depth of ALT; Sample Allelic Depth if present, otherwise based on Fillout Allelic Depth"' >> comments
                echo '# depth_sample="Read Depth (Sample)"' >> comments
                echo '# ref_count_sample="Allelic Depths of REF (Sample)"' >> comments
                echo '# alt_count_sample="Allelic Depths of ALT (Sample)"' >> comments
                echo '# is_fillout="Whether the variant was present in the original Sample (FALSE) or was generated via Fillout from related samples (TRUE)"' >> comments

                cat comments > \${output_unfiltered_filename}
                cat comments > \${output_filtered_filename}

                cat "\${fillout_tmp}" >> \${output_unfiltered_filename}
                cat "\${fillout_filtered_tmp}" >> \${output_filtered_filename}
      inputs:
        unfiltered_mafs: File[]
        filtered_mafs: File[]
        # for the unfiltered file output;
        output_unfiltered_filename:
          type: string
          default: "output.maf"
        output_filtered_filename:
          type: string
          default: "output.filtered.maf"
      outputs:
        output_file:
          type: File
          outputBinding:
            glob: $(inputs.output_unfiltered_filename)
        filtered_file:
          type: File
          outputBinding:
            glob: $(inputs.output_filtered_filename)

  convert_to_portal_format:
    run: update_cBioPortal_data.cwl
    in:
      subcommand:
        valueFrom: ${ return "maf2portal"; }
      input_file: concat_with_comments/filtered_file
      output_filename:
        valueFrom: ${ return "fillout.portal.maf"; }
    out: [ output_file ]

outputs:
  output_file:
    type: File
    outputSource: concat_with_comments/output_file
  filtered_file:
    type: File
    outputSource: concat_with_comments/filtered_file
  portal_file:
    type: File
    outputSource: convert_to_portal_format/output_file
