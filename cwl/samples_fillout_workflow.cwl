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

  create_samples_list:
    doc: creates a list of sample_ids out of the samples record array for downstream use
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

  create_clinical_samples_list:
    doc: creates a list of clinical sample_ids out of the samples record array for downstream use
    in:
      samples: samples
    out: [ clinical_sample_ids ]
    run:
      class: ExpressionTool
      inputs:
        samples: "types.yml#FilloutSample[]"
      outputs:
        clinical_sample_ids: string[]
      expression: |
        ${
          var sample_ids = [];
          for ( var i in inputs.samples ){
              if (inputs.samples[i]['sample_type'] === 'clinical' ){
                sample_ids.push(inputs.samples[i]['sample_id']);
              };
            };
          return {'clinical_sample_ids': sample_ids};
        }

  create_bam_list:
    doc: gets a list of just the bam files for downstream process
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


  # NOTE: This is important; do NOT try to do complex manipulations on maf format file, do it on vcf format instead
  maf2vcf:
    doc: converts all maf input files back to vcf for downstream processing
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

  create_fillout_targets_list:
    doc: merge all the vcf files together to create the list of target regions for GetBaseCountsMultiSample for fillout
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

  # TODO: convert this to a `scatter` step that runs per-sample in parallel, then merge the outputs otherwise we will hit the command line arg length issues eventually
  # NOTE: maybe do not do this ^^^ because currently its useful to have a multi-sample vcf output, need to investigate this more
  gbcms:
    doc: run GetBaseCountsMultiSample on all the bam files against the target regions (the merged vcf from all samples)
    in:
      sample_ids: create_samples_list/sample_ids
      bam_files: create_bam_list/bam_files
      targets_vcf: create_fillout_targets_list/merged_vcf
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

  fix_labels_and_merge_vcfs:
    doc: clean up, re-annotate, and re-header the fillout vcf, also merge against the fillout targets vcf in order to add a SRC sample source tag to tell which samples each variant originated in, and pull back in each variants original quality metrics
    in:
      fillout_vcf: gbcms/output_file
      merged_vcf: create_fillout_targets_list/merged_vcf
      merged_vcf_gz: create_fillout_targets_list/merged_vcf_gz
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


  # TODO: deprecate the need for bcftools +split by making GetBaseCountsMultiSample run individually per-sample
  split_filter_vcf:
    doc: split the multi-sample vcf into individual per-sample vcf files and apply conditional filters
    scatter: sample
    in:
      sample: samples
      fillout_vcf: fix_labels_and_merge_vcfs/fillout_sources_vcf
      clinical_sample_ids: create_clinical_samples_list/clinical_sample_ids
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
                set -eu
                # INPUTS:
                input_vcf='${ return inputs.fillout_vcf.basename ; }'
                sample_type='${ return inputs.sample["sample_type"] ; }'
                sample_id='${ return inputs.sample["sample_id"] ; }'
                clinical_sample_ids='${ return inputs.clinical_sample_ids.join(" "); }'

                clinical_samples_txt=clinical_samples.txt
                filtered_fillout_filename=filtered.vcf
                unfiltered_fillout_filename=unfiltered.vcf

                # OUTPUTS:
                # file where we will save the final filtered variants for conversion to maf
                filtered_vcf="\${sample_id}.filtered.vcf"
                # final output filename
                unfiltered_vcf="\${sample_id}.vcf"

                # make a file with the list of clinical sample ID's
                touch "\${clinical_samples_txt}"
                for i in \${clinical_sample_ids}; do echo "\${i}" >> "\${clinical_samples_txt}" ; done

                # if there is at least 1 clinical sample, then apply filters
                if [[ \$(wc -l <"\${clinical_samples_txt}") -ge 1 ]]; then
                  # filter spec:
                  # Rule 1: If a mutation is in ANY clinical sample; that is if there is a mutation in any clinical sample MAF then that mutation MUST be in the fillout.
                  # Rule 2: This is the case where NONE of the clinical samples had the mutation. It has two parts
                  # Rule 2A: If ANY of the clinical samples has a
                  # t_FL_VF >= 0.10
                  # Then DO NOT Add this mutation to the fillout
                  # Rule 2B: If ALL of the clinical samples have
                  # t_FL_VF < 0.10
                  # Then put this mutation in the fillout.
                  # add or subtract a mutation for the fillout for ALL samples. There should never be partial cases where some of the samples have fill out information and not others.

                  # AD='.' : invalid AD value means the mutation was filled out for these samples
                  # FL_VF>0.1 : the filled out variant had a frequency >0.1 (10%) for these samples
                  # @clinical_samples.txt applies filters to the samples on the FORMAT tag specified
                  # make sure to use && to apply conditions to entire variant row in vcf
                  bcftools filter -e \
                  "AD[@\${clinical_samples_txt}:*]='.' && FL_VF[@\${clinical_samples_txt}]>0.1" \
                  "\${input_vcf}" > "\${filtered_fillout_filename}"

                  cp "\${input_vcf}" "\${unfiltered_fillout_filename}"
                else
                  # if there were no clinical samples supplied then there is no need for filtering
                  cp "\${input_vcf}" "\${filtered_fillout_filename}"
                  cp "\${input_vcf}" "\${unfiltered_fillout_filename}"
                fi


                # split the multi-sample vcf into per-sample vcf files
                # dirs to hold split vcf contents
                split_dir_filtered="split_filtered"
                split_dir_unfiltered="split_unfiltered"
                # the split output vcf filepath for this sample
                split_vcf_filtered="\${split_dir_filtered}/\${sample_id}.vcf"
                split_vcf_unfiltered="\${split_dir_unfiltered}/\${sample_id}.vcf"

                mkdir -p "\${split_dir_filtered}"
                mkdir -p "\${split_dir_unfiltered}"

                bcftools +split "\${filtered_fillout_filename}" --output-type v --output "\${split_dir_filtered}"
                bcftools +split "\${unfiltered_fillout_filename}" --output-type v --output "\${split_dir_unfiltered}"

                # set the final output files
                cp "\${split_vcf_filtered}" "\${filtered_vcf}"
                cp "\${split_vcf_unfiltered}" "\${unfiltered_vcf}"

      inputs:
        sample: "types.yml#FilloutSample"
        clinical_sample_ids: string[]
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

  vcf_to_maf:
    doc: converts each sample vcf back to maf format for cBioPortal and end users
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


  concat_with_comments:
    doc: combine all the individual mafs into a single maf; add comment headers; fix some values
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
    doc: make a copy of the filtered maf that is updated for cBioPortal input format requirements
    run: update_cBioPortal_data.cwl
    in:
      subcommand:
        valueFrom: ${ return "maf2portal"; }
      input_file: concat_with_comments/filtered_file
      output_filename:
        valueFrom: ${ return "fillout.portal.maf"; }
    out: [ output_file ]

  split_uncalled_variants:
    doc: split the portal fillout file into a second file to hold uncalled variants
    run: filterUncalledMutations.cwl
    in:
      input_file: convert_to_portal_format/output_file
    out: [ called_file, uncalled_file ]

outputs:
  output_file:
    type: File
    outputSource: concat_with_comments/output_file
  filtered_file:
    type: File
    outputSource: concat_with_comments/filtered_file
  portal_file:
    type: File
    outputSource: split_uncalled_variants/called_file
  uncalled_file:
    type: File
    outputSource: split_uncalled_variants/uncalled_file
