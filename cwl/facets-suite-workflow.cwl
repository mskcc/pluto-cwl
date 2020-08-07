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
    out: [ purity_seg, hisens_seg, qc_txt, gene_level_txt, arm_level_txt, output_txt, purity_rds, hisens_rds, failed_txt, stdout_txt, stderr_txt ]

  check_run_facets:
      in:
        hisens_rds: run_facets/hisens_rds
        facets_txt: run_facets/output_txt
        failed_txt: run_facets/failed_txt
        purity_seg: run_facets/purity_seg
        hisens_seg: run_facets/hisens_seg
        qc_txt: run_facets/qc_txt
        gene_level_txt: run_facets/gene_level_txt
        arm_level_txt: run_facets/arm_level_txt
        purity_rds: run_facets/purity_rds
      out: [ hisens_rds, facets_txt, single_facets_txt, purity_seg, hisens_seg, qc_txt, gene_level_txt, arm_level_txt, purity_rds ]
      run:
          class: ExpressionTool
          id: check_facets
          inputs:
            hisens_rds: File?
            facets_txt: File?
            failed_txt: File?
            purity_seg: File?
            hisens_seg: File?
            qc_txt: File?
            gene_level_txt: File?
            arm_level_txt: File?
            purity_rds: File?
          outputs:
            hisens_rds:
              type:
                type: array
                items: ['null', File]
            facets_txt:
              type:
                type: array
                items: ['null', File]
            single_facets_txt: File?
            purity_seg: File?
            hisens_seg: File?
            qc_txt: File?
            gene_level_txt: File?
            arm_level_txt: File?
            purity_rds: File?
          expression: "${ var output_object = {};
            if(inputs.failed_txt){
              output_object['hisens_rds'] = [];
              output_object['facets_txt'] = [];
              output_object['single_facets_txt'] = null;
              output_object['purity_seg'] = null;
              output_object['hisens_seg'] = null;
              output_object['qc_txt'] = null;
              output_object['gene_level_txt'] = null;
              output_object['arm_level_txt'] = null;
              output_object['purity_rds'] = null;
            }
            else{
              if (inputs.hisens_rds && inputs.facets_txt){
                output_object['hisens_rds'] = [inputs.hisens_rds];
                output_object['facets_txt'] = [inputs.facets_txt];
                output_object['single_facets_txt'] = inputs.facets_txt;
                output_object['purity_seg'] = inputs.purity_seg;
                output_object['hisens_seg'] = inputs.hisens_seg;
                output_object['qc_txt'] = inputs.qc_txt;
                output_object['gene_level_txt'] = inputs.gene_level_txt;
                output_object['arm_level_txt'] = inputs.arm_level_txt;
                output_object['purity_rds'] = inputs.purity_rds;
              }
              else{
                output_object['hisens_rds'] = [];
                output_object['facets_txt'] = [];
                output_object['single_facets_txt'] = null;
                output_object['purity_seg'] = null;
                output_object['hisens_seg'] = null;
                output_object['qc_txt'] = null;
                output_object['gene_level_txt'] = null;
                output_object['arm_level_txt'] = null;
                output_object['purity_rds'] = null;
              }
            }
            return output_object;
          }"

  # need to run in legacy mode to get the .cncf files for downstream usages
  run_facets_legacy:
    run: run-facets-legacy-wrapper.cwl
    in:
      snp_pileup: snp_pileup
      sample_id: tumor_id
    out:
      [ hisens_cncf_txt, failed_txt, stdout_txt, stderr_txt ]

  check_run_facets_legacy:
      in:
        hisens_cncf_txt: run_facets_legacy/hisens_cncf_txt
        failed_txt: run_facets_legacy/failed_txt
      out: [ hisens_cncf_txt ]
      run:
          class: ExpressionTool
          id: check_facets_legacy
          inputs:
            hisens_cncf_txt: File?
            failed_txt: File?
          outputs:
            hisens_cncf_txt: File?
          expression: "${ var output_object = {};
            if(inputs.failed_txt){
              output_object['hisens_cncf_txt'] = null;
            }
            else{
              if (inputs.hisens_cncf_txt){
                output_object['hisens_cncf_txt'] = inputs.hisens_cncf_txt;
              }
              else{
                output_object['hisens_cncf_txt'] = null;
              }
            }
            return output_object;
          }"

  annotate_maf:
    run: annotate-maf-wrapper.cwl
    in:
      pair_id: pair_id
      maf_file: pair_maf
      facets_rds: check_run_facets/hisens_rds
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.maf
    scatter: [facets_rds]
    scatterMethod: dotproduct
    out: [output_file, failed_txt, stdout_txt, stderr_txt ]

  check_annotate_maf:
      in:
        output_file: annotate_maf/output_file
        failed_txt: annotate_maf/failed_txt
      out: [ output_file ]
      run:
          class: ExpressionTool
          id: check_annotate_maf
          inputs:
            output_file:
              type:
                type: array
                items: ['null', File]
            failed_txt:
              type:
                type: array
                items: ['null', File]
          outputs:
            output_file:
              type:
                type: array
                items: ['null', File]
          expression: "${ var output_object = {};
            if(inputs.failed_txt && inputs.failed_txt[0]){
              output_object['output_file'] = [];
            }
            else{
              if (inputs.output_file && inputs.output_file[0]){
                output_object['output_file'] = inputs.output_file;
              }
              else{
                output_object['output_file'] = [];
              }
            }
            return output_object;
          }"

  # need to apply some extra column labels to the facets suite .txt file for downstream ease of use
  label_facets_txt_tumor:
    run: paste-col.cwl
    in:
      tumor_id: tumor_id
      input_file: check_run_facets/facets_txt
      output_filename:
        valueFrom: $(inputs.tumor_id).tumor.txt
      header:
        valueFrom: ${ return "tumor"; }
      value:
        valueFrom: $(inputs.tumor_id)
    scatter: [input_file]
    scatterMethod: dotproduct
    out:
      [ output_file, failed_txt, stdout_txt, stderr_txt ]

  check_label_facets_txt_tumor:
      in:
        output_file: label_facets_txt_tumor/output_file
        failed_txt: label_facets_txt_tumor/failed_txt
      out: [ output_file ]
      run:
          class: ExpressionTool
          id: check_label_facets_txt_tumor
          inputs:
            output_file:
              type:
                type: array
                items: ['null', File]
            failed_txt:
              type:
                type: array
                items: ['null', File]
          outputs:
            output_file:
              type:
                type: array
                items: ['null', File]
          expression: "${ var output_object = {};
            if(inputs.failed_txt && inputs.failed_txt[0]){
              output_object['output_file'] = [];
            }
            else{
              if (inputs.output_file && inputs.output_file[0]){
                output_object['output_file'] = inputs.output_file;
              }
              else{
                output_object['output_file'] = [];
              }
            }
            return output_object;
          }"

  label_facets_txt_normal:
    run: paste-col.cwl
    in:
      tumor_id: tumor_id
      normal_id: normal_id
      input_file: check_label_facets_txt_tumor/output_file
      output_filename:
        valueFrom: $(inputs.tumor_id).txt
      header:
        valueFrom: ${ return "normal"; }
      value:
        valueFrom: $(inputs.normal_id)
    scatter: [input_file]
    scatterMethod: dotproduct
    out:
      [ output_file, failed_txt, stdout_txt, stderr_txt ]

  check_label_facets_txt_normal:
      in:
        output_file: label_facets_txt_normal/output_file
        failed_txt: label_facets_txt_normal/failed_txt
      out: [ output_file ]
      run:
          class: ExpressionTool
          id: check_label_facets_txt_normal
          inputs:
            output_file:
              type:
                type: array
                items: ['null', File]
            failed_txt:
              type:
                type: array
                items: ['null', File]
          outputs:
            output_file:
              type:
                type: array
                items: ['null', File]
          expression: "${ var output_object = {};
            if(inputs.failed_txt && inputs.failed_txt[0]){
              output_object['output_file'] = [];
            }
            else{
              if (inputs.output_file && inputs.output_file[0]){
                output_object['output_file'] = inputs.output_file;
              }
              else{
                output_object['output_file'] = [];
              }
            }
            return output_object;
          }"

  # need to apply some extra column labels to the maf file for downstream ease of use
  label_maf_sample:
    run: paste-col.cwl
    in:
      pair_id: pair_id
      tumor_id: tumor_id
      input_file: check_annotate_maf/output_file
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.sample.maf
      header:
        valueFrom: ${ return "sample"; }
      value:
        valueFrom: $(inputs.tumor_id)
    scatter: [input_file]
    scatterMethod: dotproduct
    out:
      [ output_file, failed_txt, stdout_txt, stderr_txt ]

  check_label_maf_sample:
      in:
        output_file: label_maf_sample/output_file
        failed_txt: label_maf_sample/failed_txt
      out: [ output_file ]
      run:
          class: ExpressionTool
          id: check_label_maf_sample
          inputs:
            output_file:
              type:
                type: array
                items: ['null', File]
            failed_txt:
              type:
                type: array
                items: ['null', File]
          outputs:
            output_file:
              type:
                type: array
                items: ['null', File]
          expression: "${ var output_object = {};
            if(inputs.failed_txt && inputs.failed_txt[0]){
              output_object['output_file'] = [];
            }
            else{
              if (inputs.output_file && inputs.output_file[0]){
                output_object['output_file'] = inputs.output_file;
              }
              else{
                output_object['output_file'] = [];
              }
            }
            return output_object;
          }"

  label_maf_normal:
    run: paste-col.cwl
    in:
      pair_id: pair_id
      normal_id: normal_id
      input_file: check_label_maf_sample/output_file
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.sample.normal.maf
      header:
        valueFrom: ${ return "normal"; }
      value:
        valueFrom: $(inputs.normal_id)
    scatter: [input_file]
    scatterMethod: dotproduct
    out:
      [ output_file, failed_txt, stdout_txt, stderr_txt ]

  check_label_maf_normal:
      in:
        output_file: label_maf_normal/output_file
        failed_txt: label_maf_normal/failed_txt
      out: [ output_file ]
      run:
          class: ExpressionTool
          id: check_label_maf_normal
          inputs:
            output_file:
              type:
                type: array
                items: ['null', File]
            failed_txt:
              type:
                type: array
                items: ['null', File]
          outputs:
            output_file:
              type:
                type: array
                items: ['null', File]
          expression: "${ var output_object = {};
            if(inputs.failed_txt && inputs.failed_txt[0]){
              output_object['output_file'] = [];
            }
            else{
              if (inputs.output_file && inputs.output_file[0]){
                output_object['output_file'] = inputs.output_file;
              }
              else{
                output_object['output_file'] = [];
              }
            }
            return output_object;
          }"

  # need to add some extra columns to the maf file from the facets output for use with cBioPortal
  update_maf:
    run: update_cBioPortal_data.cwl
    in:
      pair_id: pair_id
      subcommand:
        valueFrom: ${ return "mutations"; }
      input_file: check_label_maf_normal/output_file
      output_filename:
        valueFrom: $(inputs.pair_id)_hisens.ccf.portal.maf
      facets_txt: check_run_facets/single_facets_txt
    scatter: [input_file]
    scatterMethod: dotproduct
    out:
      [ output_file, failed_txt, stdout_txt, stderr_txt ]

  check_update_maf:
      in:
        output_file: update_maf/output_file
        failed_txt: update_maf/failed_txt
      out: [ output_file ]
      run:
          class: ExpressionTool
          id: check_update_maf
          inputs:
            output_file:
              type:
                type: array
                items: ['null', File]
            failed_txt:
              type:
                type: array
                items: ['null', File]
          outputs:
            output_file:
              type:
                type: array
                items: ['null', File]
          expression: "${ var output_object = {};
            if(inputs.failed_txt && inputs.failed_txt[0]){
              output_object['output_file'] = [];
            }
            else{
              if (inputs.output_file && inputs.output_file[0]){
                output_object['output_file'] = inputs.output_file;
              }
              else{
                output_object['output_file'] = [];
              }
            }
            return output_object;
          }"

  check_results:
      in:
        pair_id: pair_id
        hisens_cncf_txt: check_run_facets_legacy/hisens_cncf_txt
        purity_seg: check_run_facets/purity_seg
        hisens_seg: check_run_facets/hisens_seg
        qc_txt: check_run_facets/qc_txt
        gene_level_txt: check_run_facets/gene_level_txt
        arm_level_txt: check_run_facets/arm_level_txt
        facets_txt: check_label_facets_txt_normal/output_file
        purity_rds: check_run_facets/purity_rds
        hisens_rds: check_run_facets/hisens_rds
        annotated_maf: check_update_maf/output_file
        log_files:
          source: [ run_facets/stdout_txt,run_facets_legacy/stdout_txt,annotate_maf/stdout_txt,label_facets_txt_tumor/stdout_txt,label_facets_txt_normal/stdout_txt,label_maf_sample/stdout_txt,label_maf_normal/stdout_txt,update_maf/stdout_txt,run_facets/stderr_txt,run_facets_legacy/stderr_txt,annotate_maf/stderr_txt,label_facets_txt_tumor/stderr_txt,label_facets_txt_normal/stderr_txt,label_maf_sample/stderr_txt,label_maf_normal/stderr_txt,update_maf/stderr_txt]
          linkMerge: merge_flattened
      out: [ hisens_cncf_txt,purity_seg,hisens_seg,qc_txt,gene_level_txt,arm_level_txt,facets_txt,purity_rds,hisens_rds,annotated_maf,output_dir,results_passed ]
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
            hisens_rds:
              type:
                type: array
                items: ['null', File]
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
            output_dir: Directory
            results_passed: boolean
          expression: "${ var output_object = {};
            var results_passed = true;
            var facets_files = [];
            for(var key in inputs){
              var output_value = inputs[key];
              if (key == 'annotated_maf' || key == 'facets_txt' || key == 'hisens_rds'){
                if (!Array.isArray(output_value) || !output_value.length) {
                  results_passed = false;
                  output_object[key] = null;
                }
                else{
                  output_object[key] = output_value[0];
                  facets_files.push(output_value[0]);
                }
              }
              else if (key != 'log_files'){
                if ( ! output_value || Object.keys(output_value).length === 0 ){
                  results_passed = false;
                  output_object[key] = null;
                }
                else{
                  output_object[key] = output_value;
                  if ( key != 'pair_id' ){
                    facets_files.push(output_value);
                  }
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
  output_dir:
    type: Directory
    outputSource: check_results/output_dir
  results_passed:
    type: boolean
    outputSource: check_results/results_passed
