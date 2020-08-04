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
      facets_txt: check_facets/facets_txt
    scatter: [input_file,facets_txt]
    scatterMethod: dotproduct
    out:
      [ output_file, stdout_txt, stderr_txt ]

  check_results:
      in:
        pair_id: pair_id
        hisens_cncf_txt: run_facets_legacy/hisens_cncf_txt
        purity_seg: run_facets/purity_seg
        hisens_seg: run_facets/hisens_seg
        qc_txt: run_facets/qc_txt
        gene_level_txt: run_facets/gene_level_txt
        arm_level_txt: run_facets/arm_level_txt
        facets_txt: label_facets_txt_normal/output_file
        purity_rds: run_facets/purity_rds
        hisens_rds: run_facets/hisens_rds
        annotated_maf: update_maf/output_file
        log_files:
          source: [ run_facets/stdout_txt,run_facets_legacy/stdout_txt,annotate_maf/stdout_txt,label_facets_txt_tumor/stdout_txt,label_facets_txt_normal/stdout_txt,label_maf_sample/stdout_txt,label_maf_normal/stdout_txt,update_maf/stdout_txt,run_facets/stderr_txt,run_facets_legacy/stderr_txt,annotate_maf/stderr_txt,label_facets_txt_tumor/stderr_txt,label_facets_txt_normal/stderr_txt,label_maf_sample/stderr_txt,label_maf_normal/stderr_txt,update_maf/stderr_txt]
          linkMerge: merge_flattened
      out: [ hisens_cncf_txt,purity_seg,hisens_seg,qc_txt,gene_level_txt,arm_level_txt,facets_txt,purity_rds,hisens_rds,annotated_maf,log_files, results_passed ]
      run:
          class: ExpressionTool
          id: check_results
          inputs:
            pair_id: string
            hisens_cncf_txt: File?
            purity_seg: File?
            hisens_seg: File?
            qc_txt: File?
            gene_level_txt: File?
            arm_level_txt: File?
            facets_txt:
              type:
                type: array
                items: ['null', File]
            purity_rds: File?
            hisens_rds: File?
            annotated_maf:
              type:
                type: array
                items: ['null', File]
            log_files:
              type:
                type: array
                items: ['null', File]
          outputs:
            hisens_cncf_txt: File?
            purity_seg: File?
            hisens_seg: File?
            qc_txt: File?
            gene_level_txt: File?
            arm_level_txt: File?
            facets_txt: File?
            purity_rds: File?
            hisens_rds: File?
            annotated_maf: File?
            log_files: Directory?
            results_passed: boolean
          expression: "${ var output_object = {};
            var results_passed = true;
            for(var key in inputs){
              var output_value = inputs[key];
              if (key == 'annotated_maf' || key == 'facets_txt'){
                if (!Array.isArray(output_value) || !output_value.length) {
                  results_passed = false;
                  output_object[key] = null;
                }
                else{
                  output_object[key] = output_value[0];
                }
              }
              else if (key == 'log_files'){
                output_object['log_files'] = {
                  'class': 'Directory',
                  'basename': inputs.pair_id,
                  'listing': output_value
                }
              }
              else {
                if ( ! output_value || Object.keys(output_value).length === 0 ){
                  results_passed = false;
                  output_object[key] = null
                }
                else{
                  output_object[key] = output_value;
                }
              }

            }

            output_object['results_passed'] = results_passed;
            return output_object;
          }"


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
