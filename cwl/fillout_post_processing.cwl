#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
id: fillout_post_processing
label: fillout_post_processing
doc: "
Workflow to apply post-processing to fillout files after GetBaseCountsMultiSample has been applied.
This workflow can also be used on files that did not have GetBaseCountsMultiSample (fillout) applied e.g. singleton's
"

requirements:
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
  - $import: types.yml

inputs:
  samples: "types.yml#FilloutNoMafsample[]"
  fillout_vcf:
    type: File
    doc: the multi-sample .vcf file with all the mutations to be processed in the workflow, MUST contain all of the sample_ids from samples entry!
  clinical_sample_ids:
    type: string[]
    doc: a list of the sample ID's inside the fillout_vcf that are for clinical samples
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

steps:
  fillout_clinical_filter:
    doc: split the multi-sample vcf into individual per-sample vcf files and apply conditional filters
    run: fillout_clinical_filter.cwl
    scatter: sample
    in:
      sample: samples
      fillout_vcf: fillout_vcf
      clinical_sample_ids: clinical_sample_ids
    out: [ sample ] # NOTE: adds unfiltered_vcf and filtered_vcf to sample entry !!

  vcf_to_maf:
    doc: converts each sample vcf back to maf format for cBioPortal and end users
    scatter: sample
    in:
      sample: fillout_clinical_filter/sample
      ref_fasta: ref_fasta
      exac_filter: exac_filter
    out: [ unfiltered_maf, filtered_maf, sample ]
    run:
      class: CommandLineTool
      id: vcf_to_maf
      label: vcf_to_maf
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
        sample: "types.yml#FilloutNoMafsample"
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
          type: "types.yml#FilloutNoMafsample"
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
      id: concat_with_comments
      label: concat_with_comments
      baseCommand: [ "bash", "run.sh" ]
      requirements:
        DockerRequirement:
          dockerPull: mskcc/helix_filters_01:21.7.1
        InitialWorkDirRequirement:
          listing:
            - entryname: run.sh
              entry: |-
                set -eux
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
    out: [ called_file, uncalled_file ] # [ data_mutations_extended.txt , data_mutations_uncalled.txt]

outputs:
  output_file:
    type: File
    outputSource: concat_with_comments/output_file
  filtered_file:
    type: File
    outputSource: concat_with_comments/filtered_file
  portal_file:
    doc: data_mutations_extended.txt
    type: File
    outputSource: split_uncalled_variants/called_file
  uncalled_file:
    doc: data_mutations_uncalled.txt
    type: File
    outputSource: split_uncalled_variants/uncalled_file
