#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: fillout_clinical_filter
label: fillout_clinical_filter
doc: "
TODO: deprecate the need for bcftools +split by making GetBaseCountsMultiSample run individually per-sample
NOTE: maybe dont do that for ... reasons <- having single multi-sample vcf is actually easier in some ways I think

NOTE: requires AD and FL_VF FORMAT fields in the vcf!!
"

baseCommand: [ "bash", "run.sh" ]

requirements:
  - $import: types.yml
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: mskcc/helix_filters_01:21.4.1
  - class: InitialWorkDirRequirement
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
  clinical_sample_ids: 
    type: string[]
    doc: 
  fillout_vcf: 
    type: File
    doc: 

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
