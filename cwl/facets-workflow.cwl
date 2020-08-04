#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: '
Workflow for running Facets-suite on a set of tumor normal pairs

This workflow scatters over all the pairs in the input JSON to run all samples in parallel

Input JSON format
-----------------

{
  "pairs": [
    {
      "tumor_bam": {
        "class": "File",
        "path": "/test_data/bam/Tumor1.rg.md.abra.printreads.bam"
      },
      "normal_bam": {
        "class": "File",
        "path": "/test_data/bam/Normal1.rg.md.abra.printreads.bam"
      },
      "pair_maf": {
        "class": "File",
        "path": "/test_data/bam/Tumor1.Normal1.maf"
      },
      "pair_id": "Tumor1.Normal1"
    },
    {
      "tumor_bam": {
        "class": "File",
        "path": "/test_data/bam/Tumor2.rg.md.abra.printreads.bam"
      },
      "normal_bam": {
        "class": "File",
        "path": "/test_data/bam/Normal2.rg.md.abra.printreads.bam"
      },
      "pair_maf": {
        "class": "File",
        "path": "/test_data/bam/Tumor2.Normal2.maf"
      },
      "pair_id": "Tumor2.Normal2"
    }
  ]
}


Output format
-------------

output
└── facets-suite
  ├── Tumor1.Normal1.arm_level.txt
  ├── Tumor1.Normal1.gene_level.txt
  ├── Tumor1.Normal1_hisens.ccf.maf
  ├── Tumor1.Normal1_hisens.rds
  ├── Tumor1.Normal1_hisens.seg
  ├── Tumor1.Normal1_purity.rds
  ├── Tumor1.Normal1_purity.seg
  ├── Tumor1.Normal1.qc.txt
  ├── Tumor1.Normal1.snp_pileup.gz
  ├── Tumor1.Normal1.txt
  ├── Tumor2.Normal2.arm_level.txt
  ├── Tumor2.Normal2.gene_level.txt
  ├── Tumor2.Normal2_hisens.ccf.maf
  ├── Tumor2.Normal2_hisens.rds
  ├── Tumor2.Normal2_hisens.seg
  ├── Tumor2.Normal2_purity.rds
  ├── Tumor2.Normal2_purity.seg
  ├── Tumor2.Normal2.qc.txt
  ├── Tumor2.Normal2.snp_pileup.gz
  └── Tumor2.Normal2.txt
'

requirements:
  ScatterFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  pairs:
    type:
      type: array
      items:
        type: record
        fields:
          pair_maf: File
          snp_pileup: File
          pair_id: string
          tumor_id: string
          normal_id: string

steps:
  # run the facets suite wrapper set on each tumor normal pair
  facets_suite:
    run: facets-suite-workflow.cwl
    scatter: pair
    in:
      pair: pairs
      snp_pileup:
        valueFrom: ${ return inputs.pair['snp_pileup']; }
      pair_maf:
        valueFrom: ${ return inputs.pair['pair_maf']; }
      pair_id:
        valueFrom: ${ return inputs.pair['pair_id']; }
      tumor_id:
        valueFrom: ${ return inputs.pair['tumor_id']; }
      normal_id:
        valueFrom: ${ return inputs.pair['normal_id']; }
    out:
      [
      pair_id,
      purity_seg,
      hisens_seg,
      qc_txt,
      gene_level_txt,
      arm_level_txt,
      facets_txt,
      purity_rds,
      hisens_rds,
      annotated_maf,
      hisens_cncf_txt,
      log_files,
      results_passed
      ]

  check_facets_suite:
      in:
        pair_id: facets_suite/pair_id
        purity_seg: facets_suite/purity_seg
        hisens_seg: facets_suite/hisens_seg
        qc_txt: facets_suite/qc_txt
        gene_level_txt: facets_suite/gene_level_txt
        arm_level_txt: facets_suite/arm_level_txt
        facets_txt: facets_suite/facets_txt
        purity_rds: facets_suite/purity_rds
        hisens_rds: facets_suite/hisens_rds
        annotated_maf: facets_suite/annotated_maf
        hisens_cncf_txt: facets_suite/hisens_cncf_txt
        log_files: facets_suite/log_files
        results_passed: facets_suite/results_passed
      out: [ purity_seg,hisens_seg,qc_txt,gene_level_txt,arm_level_txt,facets_txt,purity_rds,hisens_rds,annotated_maf,hisens_cncf_txt,log_files,failed_pairs]
      run:
          class: ExpressionTool
          id: check_facets_suite
          inputs:
            pair_id: string[]
            purity_seg:
              type:
                type: array
                items: ['null', File]
            hisens_seg:
              type:
                type: array
                items: ['null', File]
            qc_txt:
              type:
                type: array
                items: ['null', File]
            gene_level_txt:
              type:
                type: array
                items: ['null', File]
            arm_level_txt:
              type:
                type: array
                items: ['null', File]
            facets_txt:
              type:
                type: array
                items: ['null', File]
            purity_rds:
              type:
                type: array
                items: ['null', File]
            hisens_rds:
              type:
                type: array
                items: ['null', File]
            annotated_maf:
              type:
                type: array
                items: ['null', File]
            hisens_cncf_txt:
              type:
                type: array
                items: ['null', File]
            log_files:
              type:
                type: array
                items: ['null', Directory]
            results_passed:
              type:
                type: array
                items: ['null', boolean]
          outputs:
            purity_seg:
              type:
                type: array
                items: ['null', File]
            hisens_seg:
              type:
                type: array
                items: ['null', File]
            qc_txt:
              type:
                type: array
                items: ['null', File]
            gene_level_txt:
              type:
                type: array
                items: ['null', File]
            arm_level_txt:
              type:
                type: array
                items: ['null', File]
            facets_txt:
              type:
                type: array
                items: ['null', File]
            purity_rds:
              type:
                type: array
                items: ['null', File]
            hisens_rds:
              type:
                type: array
                items: ['null', File]
            annotated_maf:
              type:
                type: array
                items: ['null', File]
            hisens_cncf_txt:
              type:
                type: array
                items: ['null', File]
            log_files: Directory?
            failed_pairs: string[]?
          expression: "${ var output_object = {};
            var passed_list = inputs.results_passed;
            var failed_list = [];
            var passed_log_dir = [];
            var failed_log_dir = [];
            output_object = {
              purity_seg: [],
              hisens_seg: [],
              qc_txt: [],
              gene_level_txt: [],
              arm_level_txt: [],
              facets_txt: [],
              purity_rds: [],
              hisens_rds: [],
              annotated_maf: [],
              hisens_cncf_txt: []
            };
            for (var pair_index in passed_list){
              if( passed_list[pair_index] == true){
                output_object['purity_seg'].push(inputs.purity_seg[pair_index]);
                output_object['hisens_seg'].push(inputs.hisens_seg[pair_index]);
                output_object['qc_txt'].push(inputs.qc_txt[pair_index]);
                output_object['gene_level_txt'].push(inputs.gene_level_txt[pair_index]);
                output_object['arm_level_txt'].push(inputs.arm_level_txt[pair_index]);
                output_object['facets_txt'].push(inputs.facets_txt[pair_index]);
                output_object['purity_rds'].push(inputs.purity_rds[pair_index]);
                output_object['hisens_rds'].push(inputs.hisens_rds[pair_index]);
                output_object['annotated_maf'].push(inputs.annotated_maf[pair_index]);
                output_object['hisens_cncf_txt'].push(inputs.hisens_cncf_txt[pair_index]);
                passed_log_dir.push(inputs.log_files[pair_index])

              }
              else{
                failed_list.push(inputs.pair_id[pair_index]);
                failed_log_dir.push(inputs.log_files[pair_index]);
              }

            }
            var success_directory = {
                  'class': 'Directory',
                  'basename': 'success',
                  'listing': passed_log_dir
            };
            var failed_directory = {
                  'class': 'Directory',
                  'basename': 'failed',
                  'listing': failed_log_dir
            };
            output_object['log_files'] = {
                  'class': 'Directory',
                  'basename': 'logs',
                  'listing': [success_directory, failed_directory]
            };
            output_object['failed_pairs'] = failed_list;
            return output_object;
          }"

outputs:
  purity_seg:
    type:
      type: array
      items: ['null', File]
    outputSource: check_facets_suite/purity_seg
  hisens_seg:
    type:
      type: array
      items: ['null', File]
    outputSource: check_facets_suite/hisens_seg
  qc_txt:
    type:
      type: array
      items: ['null', File]
    outputSource: check_facets_suite/qc_txt
  gene_level_txt:
    type:
      type: array
      items: ['null', File]
    outputSource: check_facets_suite/gene_level_txt
  arm_level_txt:
    type:
      type: array
      items: ['null', File]
    outputSource: check_facets_suite/arm_level_txt
  facets_txt:
    type:
      type: array
      items: ['null', File]
    outputSource: check_facets_suite/facets_txt
  purity_rds:
    type:
      type: array
      items: ['null', File]
    outputSource: check_facets_suite/purity_rds
  hisens_rds:
    type:
      type: array
      items: ['null', File]
    outputSource: check_facets_suite/hisens_rds
  annotated_maf:
    type:
      type: array
      items: ['null', File]
    outputSource: check_facets_suite/annotated_maf
  hisens_cncf_txt:
    type:
      type: array
      items: ['null', File]
    outputSource: check_facets_suite/hisens_cncf_txt
  log_files:
    type: Directory
    outputSource: check_facets_suite/log_files
  failed_pairs:
    type: string[]
    outputSource: check_facets_suite/failed_pairs
