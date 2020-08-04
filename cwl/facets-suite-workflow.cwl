#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: '
Workflow for running the facets suite workflow on a single tumor normal pair
'
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}

inputs:
  pair_maf: File
  snp_pileup: File
  pair_id: string
  tumor_id: string
  normal_id: string

steps:
  run_facets:
    run: run-facets-wrapper.cwl
    in:
      snp_pileup: snp_pileup
      sample_id: tumor_id
    out: [ purity_seg, hisens_seg, qc_txt, gene_level_txt, arm_level_txt, output_txt, purity_rds, hisens_rds, stdout_txt, stderr_txt ]

  # need to run in legacy mode to get the .cncf files for downstream usages
  run_facets_legacy:
    run: run-facets-legacy-wrapper.cwl
    in:
      snp_pileup: snp_pileup
      sample_id: tumor_id
    out:
      [ hisens_cncf_txt, stdout_txt, stderr_txt ]

  check_facets:
      in:
        hisens_rds: run_facets/hisens_rds
        facets_txt: run_facets/output_txt
      out: [ hisens_rds, facets_txt ]
      run:
          class: ExpressionTool
          id: check_facets
          inputs:
            hisens_rds: File?
            facets_txt: File?
          outputs:
            hisens_rds: File[]
            facets_txt: File[]
          expression: "${ var output_object = {};
            for(var key in inputs){
              var output_value = inputs[key];
              if ( !output_value ){
                output_object[key] = [];
              }
              else{
                output_object[key] = [output_value];
              }
            }
            return output_object;
          }"

  annotate_maf:
    run: annotate-maf-wrapper.cwl
    in:
      pair_id: pair_id
      maf_file: pair_maf
      facets_rds: check_facets/hisens_rds
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.maf
    scatter: [facets_rds]
    scatterMethod: dotproduct
    out: [output_file, stdout_txt, stderr_txt ]

  # need to apply some extra column labels to the facets suite .txt file for downstream ease of use
  label_facets_txt_tumor:
    run: paste-col.cwl
    in:
      tumor_id: tumor_id
      input_file: check_facets/facets_txt
      output_filename:
        valueFrom: $(inputs.tumor_id).tumor.txt
      header:
        valueFrom: ${ return "tumor"; }
      value:
        valueFrom: $(inputs.tumor_id)
    scatter: [input_file]
    scatterMethod: dotproduct
    out:
      [ output_file, stdout_txt, stderr_txt ]

  label_facets_txt_normal:
    run: paste-col.cwl
    in:
      tumor_id: tumor_id
      normal_id: normal_id
      input_file: label_facets_txt_tumor/output_file
      output_filename:
        valueFrom: $(inputs.tumor_id).txt
      header:
        valueFrom: ${ return "normal"; }
      value:
        valueFrom: $(inputs.normal_id)
    scatter: [input_file]
    scatterMethod: dotproduct
    out:
      [ output_file, stdout_txt, stderr_txt ]

  # need to apply some extra column labels to the maf file for downstream ease of use
  label_maf_sample:
    run: paste-col.cwl
    in:
      pair_id: pair_id
      tumor_id: tumor_id
      input_file: annotate_maf/output_file
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.sample.maf
      header:
        valueFrom: ${ return "sample"; }
      value:
        valueFrom: $(inputs.tumor_id)
    scatter: [input_file]
    scatterMethod: dotproduct
    out:
      [ output_file, stdout_txt, stderr_txt ]

  label_maf_normal:
    run: paste-col.cwl
    in:
      pair_id: pair_id
      normal_id: normal_id
      input_file: label_maf_sample/output_file
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.sample.normal.maf
      header:
        valueFrom: ${ return "normal"; }
      value:
        valueFrom: $(inputs.normal_id)
    scatter: [input_file]
    scatterMethod: dotproduct
    out:
      [ output_file, stdout_txt, stderr_txt ]

  # need to add some extra columns to the maf file from the facets output for use with cBioPortal
  update_maf:
    run: update_cBioPortal_data.cwl
    in:
      pair_id: pair_id
      subcommand:
        valueFrom: ${ return "mutations"; }
      input_file: label_maf_normal/output_file
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.portal.maf
      facets_txt: run_facets/output_txt
    out:
      [ output_file ]


outputs:
  pair_id:
    type: string
    outputSource: pair_id
  hisens_cncf_txt:
    type: File?
    outputSource: check_results/hisens_cncf_txt
  purity_seg:
    type: File?
    outputSource: check_results/purity_seg
  hisens_seg:
    type: File?
    outputSource: check_results/hisens_seg
  qc_txt:
    type: File?
    outputSource: check_results/qc_txt
  gene_level_txt:
    type: File?
    outputSource: check_results/gene_level_txt
  arm_level_txt:
    type: File?
    outputSource: check_results/arm_level_txt
  facets_txt:
    type: File?
    outputSource: check_results/facets_txt
  purity_rds:
    type: File?
    outputSource: check_results/purity_rds
  hisens_rds:
    type: File?
    outputSource: check_results/hisens_rds
  annotated_maf:
    type: File?
    outputSource: check_results/annotated_maf
  log_files:
    type: Directory
    outputSource: check_results/log_files
  results_passed:
    type: boolean
    outputSource: check_results/results_passed
