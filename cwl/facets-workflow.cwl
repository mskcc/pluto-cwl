#!/usr/bin/env cwl-runner

cwlVersion: v1.2
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
  ├── Tumor2.Normal2.txt
  └── logs
    ├── success
    └── failed
'

requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
  - $import: types.yml

inputs:
  pairs: "types.yml#TNMafPileupPair[]"

steps:
  # run the facets suite wrapper set on each tumor normal pair
  facets_suite:
    run: facets-suite-workflow.cwl
    doc: run the facets suite workflow separately on each TN sample pair
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
      tumor_id,
      normal_id,
      pair_id,
      purity_seg,
      hisens_seg,
      qc_txt,
      gene_level_txt,
      arm_level_txt,
      facets_txt,
      purity_rds,
      purity_png,
      hisens_rds,
      hisens_png,
      annotated_maf,
      hisens_cncf_txt,
      # output_dir,
      results_passed
      ]

  collect_pairs:
    in:
      tumor_id: facets_suite/tumor_id
      normal_id: facets_suite/normal_id
      pair_id: facets_suite/pair_id
      purity_png: facets_suite/purity_png
      purity_seg: facets_suite/purity_seg
      hisens_png: facets_suite/hisens_png
      hisens_seg: facets_suite/hisens_seg
      qc_txt: facets_suite/qc_txt
      gene_level_txt: facets_suite/gene_level_txt
      arm_level_txt: facets_suite/arm_level_txt
      facets_txt: facets_suite/facets_txt
      purity_rds: facets_suite/purity_rds
      hisens_rds: facets_suite/hisens_rds
      annotated_maf: facets_suite/annotated_maf
      hisens_cncf_txt: facets_suite/hisens_cncf_txt
    out: [ pairs ]
    run:
      class: ExpressionTool
      inputs:
        tumor_id: string[]
        normal_id: string[]
        pair_id: string[]
        purity_png:
          type:
            type: array
            items: ['null', File]
        purity_seg:
          type:
            type: array
            items: ['null', File]
        hisens_png:
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
      outputs:
        pairs:
          type: "types.yml#FacetsPair[]"
      expression: |
        ${
          var pairs = [];

          for ( var i in inputs.tumor_id ){
            var d = {
              "tumor_id": inputs.tumor_id[i],
              "normal_id": inputs.normal_id[i],
              "pair_id": inputs.pair_id[i],
              "purity_png": inputs.purity_png[i],
              "purity_seg": inputs.purity_seg[i],
              "hisens_png": inputs.hisens_png[i],
              "hisens_seg": inputs.hisens_seg[i],
              "qc_txt": inputs.qc_txt[i],
              "gene_level_txt": inputs.gene_level_txt[i],
              "arm_level_txt": inputs.arm_level_txt[i],
              "facets_txt": inputs.facets_txt[i],
              "purity_rds": inputs.purity_rds[i],
              "hisens_rds": inputs.hisens_rds[i],
              "annotated_maf": inputs.annotated_maf[i],
              "hisens_cncf_txt": inputs.hisens_cncf_txt[i]
            };
            pairs.push(d);
          };

          return {"pairs": pairs};
        }













  check_facets_suite:
      in:
        pair_id: facets_suite/pair_id
        purity_png: facets_suite/purity_png
        purity_seg: facets_suite/purity_seg
        hisens_png: facets_suite/hisens_png
        hisens_seg: facets_suite/hisens_seg
        qc_txt: facets_suite/qc_txt
        gene_level_txt: facets_suite/gene_level_txt
        arm_level_txt: facets_suite/arm_level_txt
        facets_txt: facets_suite/facets_txt
        purity_rds: facets_suite/purity_rds
        hisens_rds: facets_suite/hisens_rds
        annotated_maf: facets_suite/annotated_maf
        hisens_cncf_txt: facets_suite/hisens_cncf_txt
        # output_dir: facets_suite/output_dir
        results_passed: facets_suite/results_passed
      out: [
        purity_png,
        purity_seg,
        hisens_png,
        hisens_seg,
        qc_txt,
        gene_level_txt,
        arm_level_txt,
        facets_txt,
        purity_rds,
        hisens_rds,
        annotated_maf,
        hisens_cncf_txt,
        # output_dir,
        failed_pairs
      ]
      run:
          class: ExpressionTool
          id: check_facets_suite
          inputs:
            pair_id: string[]
            purity_png:
              type:
                type: array
                items: ['null', File]
            purity_seg:
              type:
                type: array
                items: ['null', File]
            hisens_png:
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
            # output_dir:
            #   type:
            #     type: array
            #     items: [Directory]
            results_passed:
              type:
                type: array
                items: ['null', boolean]
          outputs:
            purity_png:
              type:
                type: array
                items: ['null', File]
            purity_seg:
              type:
                type: array
                items: ['null', File]
            hisens_png:
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
            # output_dir: Directory
            failed_pairs: string[]?
          expression: "${ var output_object = {};
            var passed_list = inputs.results_passed;
            var failed_list = [];
            var output_list = [];
            output_object = {
              purity_png: [],
              purity_seg: [],
              hisens_seg: [],
              hisens_png: [],
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
                output_object['purity_png'].push(inputs.purity_png[pair_index]);
                output_object['purity_seg'].push(inputs.purity_seg[pair_index]);
                output_object['hisens_png'].push(inputs.hisens_png[pair_index]);
                output_object['hisens_seg'].push(inputs.hisens_seg[pair_index]);
                output_object['qc_txt'].push(inputs.qc_txt[pair_index]);
                output_object['gene_level_txt'].push(inputs.gene_level_txt[pair_index]);
                output_object['arm_level_txt'].push(inputs.arm_level_txt[pair_index]);
                output_object['facets_txt'].push(inputs.facets_txt[pair_index]);
                output_object['purity_rds'].push(inputs.purity_rds[pair_index]);
                output_object['hisens_rds'].push(inputs.hisens_rds[pair_index]);
                output_object['annotated_maf'].push(inputs.annotated_maf[pair_index]);
                output_object['hisens_cncf_txt'].push(inputs.hisens_cncf_txt[pair_index]);

              }
              else{
                failed_list.push(inputs.pair_id[pair_index]);
              }
              // output_list.push(inputs.output_dir[pair_index])

            }

            output_object['failed_pairs'] = failed_list;
            return output_object;
          }"
        # // output_object['output_dir'] = {
        # //      'class': 'Directory',
        # //      'basename': 'facets',
        # //      'listing': output_list
        # //};

outputs:
  pairs:
    type: "types.yml#FacetsPair[]"
    outputSource: collect_pairs/pairs
