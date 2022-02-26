#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
doc: '
Workflow for running the facets suite workflow on a single tumor normal pair

Includes handling of errors in case execution fails for the sample pair
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
    out: [
    purity_seg,
    hisens_seg,
    qc_txt,
    gene_level_txt,
    arm_level_txt,
    output_txt,
    purity_rds,
    hisens_rds,
    failed,
    stdout_txt,
    stderr_txt
    ]

  # need to run in legacy mode to get the .cncf files for downstream usages
  run_facets_legacy:
    run: run-facets-legacy-wrapper.cwl
    in:
      snp_pileup: snp_pileup
      sample_id: tumor_id
    out:
      [ hisens_cncf_txt, failed, stdout_txt, stderr_txt ]


  annotate_maf:
    run: annotate-maf-wrapper.cwl
    in:
      pair_id: pair_id
      facets_failed: run_facets/failed
      maf_file: pair_maf
      facets_rds: run_facets/hisens_rds
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.maf
    when: $(inputs.facets_failed == false && inputs.facets_rds != null)
    out: [output_file, failed, stdout_txt, stderr_txt ]

  # need to apply some extra column labels to the facets suite .txt file for downstream ease of use
  label_facets_txt_tumor:
    run: paste-col-wrapper.cwl
    in:
      tumor_id: tumor_id
      facets_failed: run_facets/failed
      input_file: run_facets/output_txt
      output_filename:
        valueFrom: $(inputs.tumor_id).tumor.txt
      header:
        valueFrom: ${ return "tumor"; }
      value:
        valueFrom: $(inputs.tumor_id)
    when: $(inputs.facets_failed == false && inputs.input_file != null )
    out:
      [ output_file, failed, stdout_txt, stderr_txt ]

  label_facets_txt_normal:
    run: paste-col-wrapper.cwl
    in:
      tumor_id: tumor_id
      normal_id: normal_id
      label_facets_failed: label_facets_txt_tumor/failed
      input_file: label_facets_txt_tumor/output_file
      output_filename:
        valueFrom: $(inputs.tumor_id).txt
      header:
        valueFrom: ${ return "normal"; }
      value:
        valueFrom: $(inputs.normal_id)
    when: $(inputs.label_facets_failed == false && inputs.input_file != null )
    out:
      [ output_file, failed, stdout_txt, stderr_txt ]

  # need to apply some extra column labels to the maf file for downstream ease of use
  label_maf_sample:
    run: paste-col-wrapper.cwl
    in:
      pair_id: pair_id
      tumor_id: tumor_id
      annotate_maf_failed: annotate_maf/failed
      input_file: annotate_maf/output_file
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.sample.maf
      header:
        valueFrom: ${ return "sample"; }
      value:
        valueFrom: $(inputs.tumor_id)
    when: $(inputs.annotate_maf_failed == false && inputs.input_file != null )
    out:
      [ output_file, failed, stdout_txt, stderr_txt ]  

  label_maf_normal:
    run: paste-col-wrapper.cwl
    in:
      pair_id: pair_id
      normal_id: normal_id
      label_maf_sample_failed: label_maf_sample/failed
      input_file: label_maf_sample/output_file
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.sample.normal.maf
      header:
        valueFrom: ${ return "normal"; }
      value:
        valueFrom: $(inputs.normal_id)
    when: $(inputs.label_maf_sample_failed == false && inputs.input_file != null )
    out:
      [ output_file, failed, stdout_txt, stderr_txt ]

  # need to add some extra columns to the maf file from the facets output for use with cBioPortal
  update_maf:
    run: update_cBioPortal_data.cwl
    in:
      pair_id: pair_id
      subcommand:
        valueFrom: ${ return "mutations"; }
      label_maf_normal_failed: label_maf_normal/failed
      input_file: label_maf_normal/output_file
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.portal.maf
      facets_txt: run_facets/output_txt
      facets_failed: run_facets/failed
    when: $( inputs.label_maf_normal_failed == false && inputs.facets_failed == false && inputs.input_file != null )
    out:
      [ output_file, failed, stdout_txt, stderr_txt ]

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
          source: [
            run_facets/stdout_txt,
            run_facets_legacy/stdout_txt,
            annotate_maf/stdout_txt,
            label_facets_txt_tumor/stdout_txt,
            label_facets_txt_normal/stdout_txt,
            label_maf_sample/stdout_txt,
            label_maf_normal/stdout_txt,
            update_maf/stdout_txt,
            run_facets/stderr_txt,
            run_facets_legacy/stderr_txt,
            annotate_maf/stderr_txt,
            label_facets_txt_tumor/stderr_txt,
            label_facets_txt_normal/stderr_txt,
            label_maf_sample/stderr_txt,
            label_maf_normal/stderr_txt,
            update_maf/stderr_txt
            ]
          pickValue: all_non_null
        results:
          source: [
            run_facets/failed,
            run_facets_legacy/failed,
            annotate_maf/failed,
            label_facets_txt_tumor/failed,
            label_facets_txt_normal/failed,
            label_maf_sample/failed,
            label_maf_normal/failed,
            update_maf/failed
          ]
          pickValue: all_non_null
      out: [
        hisens_cncf_txt,
        purity_seg,hisens_seg,
        qc_txt,
        gene_level_txt,
        arm_level_txt,
        facets_txt,
        purity_rds,
        hisens_rds,
        annotated_maf,
        output_dir,
        results_passed
        ]
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
            facets_txt: File?
            purity_rds: File?
            hisens_rds: File?
            annotated_maf: File?
            log_files:
              type:
                type: array
                items: ['null', File]
            results:
              type:
                type: array
                items: ['null', boolean]
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
            output_dir: Directory
            results_passed: boolean
          expression: "${ var output_object = {};
            var results_passed = !inputs.results.includes(true);
            if( inputs.results.length != 8 ){
              results_passed = false
            }
            var facets_files = [];
            for(var key in inputs){
              var output_value = inputs[key];
              if (key != 'log_files' || key != 'results' || key != 'pair_id'){
                if ( ! output_value ){
                  results_passed = false;
                  output_object[key] = null;
                }
                else{
                  output_object[key] = output_value;
                  facets_files.push(output_value);
                }
              }

            }

            if( results_passed == true ){
              output_object['output_dir'] = {
                  'class': 'Directory',
                  'basename': inputs.pair_id,
                  'listing': facets_files
                }

            }
            else{
              for(var key in output_object){
                output_object[key] = null;
              }
              output_object['output_dir'] = {
                  'class': 'Directory',
                  'basename': inputs.pair_id,
                  'listing': inputs.log_files
                }

            }

            output_object['results_passed'] = results_passed;
            return output_object;
          }"


outputs:
  pair_id:
    type: string
    outputSource: pair_id
  hisens_cncf_txt: # Tumor1.Normal1_hisens.cncf.txt ; from legacy facets output
    type: File?
    outputSource: check_results/hisens_cncf_txt
  purity_seg: # Tumor1.Normal1_purity.seg
    type: File?
    outputSource: check_results/purity_seg
  hisens_seg: # Tumor1.Normal1_hisens.seg
    type: File?
    outputSource: check_results/hisens_seg
  qc_txt: # Tumor1.Normal1.qc.txt
    type: File?
    outputSource: check_results/qc_txt
  gene_level_txt: # Tumor1.Normal1.gene_level.txt
    type: File?
    outputSource: check_results/gene_level_txt
  arm_level_txt: # Tumor2.Normal2.arm_level.txt
    type: File?
    outputSource: check_results/arm_level_txt
  facets_txt: # Tumor1.Normal1.txt
    type: File?
    outputSource: check_results/facets_txt
  purity_rds: # Tumor1.Normal1_purity.rds
    type: File?
    outputSource: check_results/purity_rds
  hisens_rds: # Tumor1.Normal1_hisens.rds
    type: File?
    outputSource: check_results/hisens_rds
  annotated_maf: # Tumor1.Normal1_hisens.ccf.maf ; _hisens.ccf.portal.maf
    type: File?
    outputSource: check_results/annotated_maf
  output_dir:
    type: Directory
    outputSource: check_results/output_dir
  results_passed:
    type: boolean
    outputSource: check_results/results_passed
