#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
doc: "
CWL workflow for generating Roslin / Argos post pipeline analysis files and cBioPortal data and metadata files

This workflow includes Facets and Facets Suite usages

Inputs
------

The following parameters are required:

project_id
project_pi
request_pi
project_short_name
project_name
project_description
cancer_type
cancer_study_identifier
argos_version_string
helix_filter_version
is_impact
extra_pi_groups
pairs


The following filenames are required:

analysis_mutations_filename
analysis_gene_cna_filename
analysis_sv_filename
analysis_segment_cna_filename
cbio_segment_data_filename
cbio_meta_cna_segments_filename

The following filenames have default values and are optional:

cbio_mutation_data_filename
cbio_cna_data_filename
cbio_fusion_data_filename
cbio_clinical_patient_data_filename
cbio_clinical_sample_data_filename
cbio_clinical_sample_meta_filename
cbio_clinical_patient_meta_filename
cbio_meta_study_filename
cbio_meta_cna_filename
cbio_meta_fusions_filename
cbio_meta_mutations_filename
cbio_cases_all_filename
cbio_cases_cnaseq_filename
cbio_cases_cna_filename
cbio_cases_sequenced_filename

Output
------

Workflow output should look like this:

output
├── analysis
│   ├── <project_id>.gene.cna.txt
│   ├── <project_id>.muts.maf
│   ├── <project_id>.seg.cna.txt
│   └── <project_id>.svs.maf
├── facets
│   ├── <tumor_id>.<normal_id> (passed)
│   │   └── <facets_files>
│   └── <tumor_id>.<normal_id> (failed)
│       └── <log_files>
└── portal
    ├── case_list
    │   ├── cases_all.txt
    │   ├── cases_cnaseq.txt
    │   ├── cases_cna.txt
    │   └── cases_sequenced.txt
    ├── data_clinical_patient.txt
    ├── data_clinical_sample.txt
    ├── data_CNA.ascna.txt
    ├── data_CNA.scna.txt
    ├── data_CNA.txt
    ├── data_sv.txt
    ├── data_fusions.txt
    ├── data_mutations_extended.txt
    ├── meta_clinical_patient.txt
    ├── meta_clinical_sample.txt
    ├── meta_CNA.txt
    ├── meta_sv.txt
    ├── meta_fusions.txt
    ├── meta_mutations_extended.txt
    ├── meta_study.txt
    ├── <project_id>_data_cna_hg19.seg
    └── <project_id>_meta_cna_hg19_seg.txt
"

requirements:
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
  - $import: types.yml


inputs:
  pairs:
    doc: list of tumor normal sample pairs
    type: "types.yml#TNMafPileupPair[]"
  project_id:
    type: string
    doc: "unique identifier for the project"
  project_pi:
    type: string
    doc: "principle investigator for the project"
  request_pi:
    type: string
    doc: "principle investigator who requested the project"
  project_short_name:
    type: string
    doc: "a short name for the project in cBioPortal"
  project_name:
    type: string
    doc: "a formal name for the project"
  project_description:
    type: string
    doc: "a description of the project"
  cancer_type:
    type: string
    doc: "the type of cancer used in the project"
  cancer_study_identifier:
    type: string
    doc: "a study identifier for the project to use in cBioPortal"
  argos_version_string:
    type: string
    doc: "the version label of Roslin / Argos used to run the project analysis"
  helix_filter_version:
    type: string
    doc: "the version label of this helix filter repo (git describe --all --long)"
  is_impact:
    default: True
    type: boolean
    doc: "whether or not the project is an IMPACT project; should be the value 'True' if so, otherwise any other value means 'False'"
  # TODO: this shouild actually be type: string[]
  extra_pi_groups:
    type: [ "null", string]
    default: null
    doc: "a list of other groups to be associated with the project in cBioPortal"
  analysis_segment_cna_filename:
    type: string
    doc: "<project_id>.seg.cna.txt"
  analysis_sv_filename:
    type: string
    doc: "<project_id>.svs.maf"
  analysis_gene_cna_filename:
    type: string
    doc: "<project_id>.gene.cna.txt"
  analysis_mutations_filename:
    type: string
    doc: "<project_id>.muts.maf"
  analysis_mutations_share_filename:
    type: string
    doc: "<project_id>.muts.share.maf"
  cbio_segment_data_filename:
    type: string
    doc: "<project_id>_data_cna_hg19.seg"
  cbio_meta_cna_segments_filename:
    type: string
    doc: "<project_id>_meta_cna_hg19_seg.txt"
  cbio_cases_sequenced_filename:
    type: string
    doc: ""
    default: cases_sequenced.txt
  cbio_cases_cna_filename:
    type: string
    default: cases_cna.txt
    doc: ""
  cbio_cases_cnaseq_filename:
    type: string
    default: cases_cnaseq.txt
    doc: ""
  cbio_cases_all_filename:
    type: string
    default: cases_all.txt
    doc: ""
  cbio_meta_mutations_filename:
    type: string
    default: meta_mutations_extended.txt
    doc: ""
  cbio_meta_fusions_filename:
    type: string
    default: meta_fusions.txt
    doc: ""
  cbio_meta_cna_filename:
    type: string
    default: meta_CNA.txt
    doc: ""
  cbio_meta_study_filename:
    type: string
    default: meta_study.txt
    doc: ""
  cbio_clinical_patient_meta_filename:
    type: string
    default: meta_clinical_patient.txt
    doc: ""
  cbio_clinical_sample_meta_filename:
    type: string
    default: meta_clinical_sample.txt
    doc: ""
  cbio_clinical_sample_data_filename:
    type: string
    default: data_clinical_sample.txt
    doc: ""
  cbio_clinical_patient_data_filename:
    type: string
    default: data_clinical_patient.txt
    doc: ""
  cbio_fusion_data_filename:
    type: string
    default: data_fusions.txt
    doc: ""
  cbio_mutation_data_filename:
    type: string
    default: data_mutations_extended.txt
    doc: ""
  cbio_cna_data_filename:
    type: string
    default: data_CNA.txt
    doc: ""
  cbio_cna_ascna_data_filename:
    type: string
    default: data_CNA.ascna.txt
    doc: ""
  cbio_cna_scna_data_filename:
    type: string
    default: data_CNA.scna.txt
    doc: ""
  cbio_sv_data_filename:
    type: string
    default: data_sv.txt
    doc: ""
  cbio_meta_sv_filename:
    type: string
    default: meta_sv.txt
    doc: ""
  mutation_svs_txt_files:
    type: File[]
    doc: "*.svs.pass.vep.portal.txt"
  mutation_svs_maf_files:
    type: File[]
    doc: "*.svs.pass.vep.maf"
  targets_list:
    type: File
  known_fusions_file:
    type: File
  data_clinical_file:
    type: File
  sample_summary_file:
    type:
    - "null"
    - File
  IMPACT_gene_list:
    type: File
    doc: "TSV file with gene labels and corresponding impact assays"
  assay_coverage:
    type: string
    doc: "genome_coverage value; amount of the genome in bp covered by the assay"
  microsatellites_file:
    type: File
    doc: "Microsatellites list file to use with MSI Sensor"
  # NOTE: these two arrays of File with secondaryFiles should eventually be merged directly into the `pairs` record array
  # after upgrading Toil to support cwlVersion 1.1
  normal_bam_files:
    type:
        type: array
        items: File
    doc: "Array of normal bam files. Must match the same order of sample pairs in 'pairs' input field"
    secondaryFiles:
        - ^.bai
  tumor_bam_files:
    type:
        type: array
        items: File
    doc: "Array of tumor bam files. Must match the same order of sample pairs in 'pairs' input field"
    secondaryFiles:
        - ^.bai
  extra_cna_files:
    doc: "Extra CNA data files to be merged in with the portal CNA data"
    type:
    - "null"
    - File[]
  extra_sample_ids:
    doc: Extra sample ids that should be included in case list files
    type:
    - string[]
    - "null"


steps:
  run_facets:
    doc: run the Facets Suite workflow
    run: facets-workflow.cwl
    in:
      pairs: pairs
    out: [ pairs ]

  gather_facets_files:
    doc: gather some files from facets output pairs in order to pass to downstream steps as file lists
    in:
      pairs: run_facets/pairs
    out: [ annotated_mafs, facets_txts, hisens_cncf_txts, hisens_segs ]
    run:
      class: ExpressionTool
      inputs:
        pairs:
          type: "types.yml#FacetsPair[]"
      outputs:
        annotated_mafs: File[]
        facets_txts: File[]
        hisens_cncf_txts: File[]
        hisens_segs: File[]
      expression: |
        ${
          var annotated_mafs = [];
          var facets_txts = [];
          var hisens_cncf_txts = [];
          var hisens_segs = [];

          for ( var i in inputs.pairs ){
            annotated_mafs.push(inputs.pairs[i].annotated_maf)
            facets_txts.push(inputs.pairs[i].facets_txt)
            hisens_cncf_txts.push(inputs.pairs[i].hisens_cncf_txt)
            hisens_segs.push(inputs.pairs[i].hisens_seg)
          };

          return {
            "annotated_mafs": annotated_mafs,
            "facets_txts": facets_txts,
            "hisens_cncf_txts": hisens_cncf_txts,
            "hisens_segs": hisens_segs,
          };
        }


  concat_facets_maf:
    doc: make a concatenated maf file for merging with portal maf
    run: concat-tables.cwl
    in:
      input_files: gather_facets_files/annotated_mafs
      output_filename:
        valueFrom: ${ return "facets.maf"; }
    out:
      [ output_file ]



  run_analysis_workflow:
    doc: generate the analysis output files
    run: analysis-workflow.cwl
    in:
      analysis_segment_cna_filename: analysis_segment_cna_filename
      analysis_sv_filename: analysis_sv_filename
      analysis_gene_cna_filename: analysis_gene_cna_filename
      analysis_mutations_filename: analysis_mutations_filename
      analysis_mutations_share_filename: analysis_mutations_share_filename
      mutation_maf_files: gather_facets_files/annotated_mafs
      facets_hisens_seg_files: gather_facets_files/hisens_segs
      facets_hisens_cncf_files: gather_facets_files/hisens_cncf_txts
      mutation_svs_maf_files: mutation_svs_maf_files
      targets_list: targets_list
      argos_version_string: argos_version_string
      is_impact: is_impact
      helix_filter_version: helix_filter_version
      IMPACT_gene_list: IMPACT_gene_list
    out:
      [ analysis_dir ]



  convert_tmb_pairs:
    doc: convert the workflow input pairs record schema to the TMB analysis pair schema
    run: convert_TNMafPileupPair_to_TMBInputPair.cwl
    in:
      pairs: pairs
    out: [ pairs ]

  run_tmb_workflow:
    doc: run the TMB workflow
    run: tmb_workflow.cwl
    in:
      pairs: convert_tmb_pairs/pairs
      assay_coverage: assay_coverage
    out: [ pairs ]

  gather_tmb_tsvs:
    doc: create a list of just TMB tsv files for downstream processing
    in:
      pairs: run_tmb_workflow/pairs
    out: [ tmb_tsvs ]
    run:
      class: ExpressionTool
      inputs:
        pairs:
          type: "types.yml#TMBOutputPair[]"
      outputs:
        tmb_tsvs: File[]
      expression: |
        ${
          var tmb_tsvs = [];
          for ( var i in inputs.pairs ){
            tmb_tsvs.push(inputs.pairs[i].tmb_tsv);
          };
          return {"tmb_tsvs": tmb_tsvs};
        }


  convert_msi_pairs:
    run: convert_TNMafPileupPair_to_MSIInputPair.cwl
    in:
      pairs: pairs
    out: [ pairs ]

  run_msi:
    doc: run the MSI workflow
    run: msi_workflow.cwl
    in:
      microsatellites_file: microsatellites_file
      pairs: convert_msi_pairs/pairs
      normal_bam_files: normal_bam_files
      tumor_bam_files: tumor_bam_files
    out: [ pairs ]

  gather_msi_tsvs:
    doc: create list of just MSI tsv files for downstream processing
    in:
      pairs: run_msi/pairs
    out: [ msi_tsvs ]
    run:
      class: ExpressionTool
      inputs:
        pairs:
          type: "types.yml#MSIOutputPair[]"
      outputs:
        msi_tsvs: File[]
      expression: |
        ${
          var msi_tsvs = [];
          for ( var i in inputs.pairs ){
            msi_tsvs.push(inputs.pairs[i].msi_tsv);
          };
          return {"msi_tsvs": msi_tsvs};
        }











  run_portal_workflow:
    doc: generate the cBioPortal output files
    run: portal-workflow.cwl
    in:
      project_id: project_id
      project_pi: project_pi
      request_pi: request_pi
      project_short_name: project_short_name
      project_name: project_name
      project_description: project_description
      cancer_type: cancer_type
      cancer_study_identifier: cancer_study_identifier
      argos_version_string: argos_version_string
      helix_filter_version: helix_filter_version
      is_impact: is_impact
      extra_pi_groups: extra_pi_groups
      cbio_segment_data_filename: cbio_segment_data_filename
      cbio_meta_cna_segments_filename: cbio_meta_cna_segments_filename
      cbio_cases_sequenced_filename: cbio_cases_sequenced_filename
      cbio_cases_cna_filename: cbio_cases_cna_filename
      cbio_cases_cnaseq_filename: cbio_cases_cnaseq_filename
      cbio_cases_all_filename: cbio_cases_all_filename
      cbio_meta_mutations_filename: cbio_meta_mutations_filename
      cbio_meta_fusions_filename: cbio_meta_fusions_filename
      cbio_meta_cna_filename: cbio_meta_cna_filename
      cbio_meta_study_filename: cbio_meta_study_filename
      cbio_clinical_patient_meta_filename: cbio_clinical_patient_meta_filename
      cbio_clinical_sample_meta_filename: cbio_clinical_sample_meta_filename
      cbio_clinical_sample_data_filename: cbio_clinical_sample_data_filename
      cbio_clinical_patient_data_filename: cbio_clinical_patient_data_filename
      cbio_fusion_data_filename: cbio_fusion_data_filename
      cbio_mutation_data_filename: cbio_mutation_data_filename
      cbio_cna_data_filename: cbio_cna_data_filename
      cbio_cna_ascna_data_filename: cbio_cna_ascna_data_filename
      cbio_cna_scna_data_filename: cbio_cna_scna_data_filename
      mutation_maf_files: gather_facets_files/annotated_mafs
      facets_hisens_seg_files: gather_facets_files/hisens_segs
      facets_hisens_cncf_files: gather_facets_files/hisens_cncf_txts
      mutation_svs_txt_files: mutation_svs_txt_files
      targets_list: targets_list
      known_fusions_file: known_fusions_file
      data_clinical_file: data_clinical_file
      sample_summary_file: sample_summary_file
      msi_files: gather_msi_tsvs/msi_tsvs
      tmb_files: gather_tmb_tsvs/tmb_tsvs
      facets_suite_txt_files: gather_facets_files/facets_txts
      extra_sample_ids: extra_sample_ids
      extra_cna_files: extra_cna_files
      cbio_sv_data_filename: cbio_sv_data_filename
      cbio_meta_sv_filename: cbio_meta_sv_filename
    out:
      [
      portal_meta_clinical_sample_file, # meta_clinical_sample.txt
      portal_data_clinical_patient_file, # data_clinical_patient.txt
      portal_data_clinical_sample_file, # data_clinical_sample.txt
      portal_meta_study_file, # meta_study.txt
      portal_clinical_patient_meta_file, # meta_clinical_patient.txt
      portal_meta_cna_file, # meta_CNA.txt
      portal_meta_fusions_file, # meta_fusions.txt
      portal_meta_mutations_extended_file, # meta_mutations_extended.txt
      portal_meta_cna_segments_file, # <project_id>_meta_cna_hg19_seg.txt
      portal_meta_sv_file, # meta_SV.txt
      portal_cna_data_file, # data_CNA.txt
      portal_cna_ascna_file, # data_CNA.ascna.txt
      portal_muts_file, # data_mutations_extended.txt
      portal_hisens_segs, # <project_id>_data_cna_hg19.seg
      portal_fusions_data_file, # data_fusions.txt
      portal_sv_data_file, # data_SV.txt
      portal_case_list_dir,
      portal_report
      ]


  merge_maf:
    doc: need to merge the portal mutations maf with the Facets maf to get some extra information in the output
    run: update_cBioPortal_data.cwl
    in:
      subcommand:
        valueFrom: ${ return "merge_mafs"; }
      input_file: run_portal_workflow/portal_muts_file
      output_filename: cbio_mutation_data_filename
      facets_maf: concat_facets_maf/output_file
    out:
      [ output_file, failed, stdout_txt, stderr_txt ]












  # OUTPUT DIR CREATION STEPS
  # TODO: NEED TO MOVE THESE OUT OF THIS WORKFLOW INTO A DEDICATED PROCESS

  # create the "portal" directory in the output dir and put cBioPortal files in it
  make_portal_dir:
    run: put_in_dir.cwl
    in:
      portal_meta_clinical_sample_file: run_portal_workflow/portal_meta_clinical_sample_file # meta_clinical_sample.txt
      portal_data_clinical_patient_file: run_portal_workflow/portal_data_clinical_patient_file # data_clinical_patient.txt
      portal_data_clinical_sample_file: run_portal_workflow/portal_data_clinical_sample_file # data_clinical_sample.txt
      portal_meta_study_file: run_portal_workflow/portal_meta_study_file # meta_study.txt
      portal_clinical_patient_meta_file: run_portal_workflow/portal_clinical_patient_meta_file # meta_clinical_patient.txt
      portal_meta_cna_file: run_portal_workflow/portal_meta_cna_file # meta_CNA.txt
      portal_meta_fusions_file: run_portal_workflow/portal_meta_fusions_file # meta_fusions.txt
      portal_meta_mutations_extended_file: run_portal_workflow/portal_meta_mutations_extended_file # meta_mutations_extended.txt
      portal_meta_cna_segments_file: run_portal_workflow/portal_meta_cna_segments_file  # <project_id>_meta_cna_hg19_seg.txt
      portal_meta_sv_file: run_portal_workflow/portal_meta_sv_file # meta_SV.txt
      portal_cna_data_file: run_portal_workflow/portal_cna_data_file # data_CNA.txt
      portal_sv_data_file: run_portal_workflow/portal_sv_data_file # data_SV.txt
      portal_cna_ascna_file: run_portal_workflow/portal_cna_ascna_file # data_CNA.ascna.txt
      portal_muts_file: merge_maf/output_file # data_mutations_extended.txt
      portal_hisens_segs: run_portal_workflow/portal_hisens_segs # # <project_id>_data_cna_hg19.seg
      portal_fusions_data_file: run_portal_workflow/portal_fusions_data_file # data_fusions.txt
      portal_case_list_dir: run_portal_workflow/portal_case_list_dir
      portal_report: run_portal_workflow/portal_report
      output_directory_name:
        valueFrom: ${ return "portal"; }
      files:
        valueFrom: ${return [
          inputs.portal_meta_clinical_sample_file,
          inputs.portal_data_clinical_patient_file,
          inputs.portal_data_clinical_sample_file,
          inputs.portal_meta_study_file,
          inputs.portal_clinical_patient_meta_file,
          inputs.portal_meta_cna_file,
          inputs.portal_meta_fusions_file,
          inputs.portal_meta_mutations_extended_file,
          inputs.portal_meta_cna_segments_file,
          inputs.portal_meta_sv_file,
          inputs.portal_sv_data_file,
          inputs.portal_cna_data_file,
          inputs.portal_cna_ascna_file,
          inputs.portal_muts_file,
          inputs.portal_hisens_segs,
          inputs.portal_fusions_data_file,
          inputs.portal_case_list_dir,
          inputs.portal_report
          ]}
    out: [ directory ]

  make_msi_dir:
    doc:
    run: put_DirFileList_in_dir.cwl
    in:
      msi_tsvs: gather_msi_tsvs/msi_tsvs
      output_directory_name:
        valueFrom: ${ return "msi"; }
      files:
        valueFrom: ${return [
          inputs.msi_tsvs
          ]}
    out: [ directory ]

  make_tmb_dir:
    doc:
    run: put_DirFileList_in_dir.cwl
    in:
      tmb_tsvs: gather_tmb_tsvs/tmb_tsvs
      # tmb_mafs:
      output_directory_name:
        valueFrom: ${ return "tmb"; }
      files:
        valueFrom: ${return [
          inputs.tmb_tsvs
          ]}
    out: [ directory ]

  make_facets_dir:
    doc: make a single directory to hold the results for all Facets sample pair
    in:
      pairs: run_facets/pairs
    out: [ facets_dir ]
    run:
      class: Workflow
      inputs:
        pairs:
          type: "types.yml#FacetsPair[]"
      outputs:
        facets_dir:
          type: Directory
          outputSource: make_facets_dir/directory
      steps:
        make_facets_pair_dirs:
          doc: make a subdir to hold the results for one Facets sample pair
          run: put_in_dir.cwl
          scatter: pair
          in:
            pair: pairs
            output_directory_name:
              valueFrom: ${ return inputs.pair.pair_id; }
            files:
              valueFrom: ${
                  return [
                    inputs.pair.annotated_maf,
                    inputs.pair.arm_level_txt,
                    inputs.pair.facets_txt,
                    inputs.pair.gene_level_txt,
                    inputs.pair.hisens_cncf_txt,
                    inputs.pair.hisens_rds,
                    inputs.pair.hisens_seg,
                    inputs.pair.hisens_png,
                    inputs.pair.purity_rds,
                    inputs.pair.purity_seg,
                    inputs.pair.purity_png,
                    inputs.pair.qc_txt
                  ];
                }
          out: [ directory ]

        make_facets_dir:
          run: put_in_dir.cwl
          in:
            sample_dirs: make_facets_pair_dirs/directory
            output_directory_name:
              valueFrom: ${ return "facets"; }
            files:
              valueFrom: ${
                return inputs.sample_dirs;
                }
          out: [ directory ]



outputs:
  portal_dir:
    type: Directory
    outputSource: make_portal_dir/directory

  analysis_dir:
    type: Directory
    outputSource: run_analysis_workflow/analysis_dir

  facets_dir:
    type: Directory
    outputSource: make_facets_dir/facets_dir

  msi_dir:
    type: Directory
    outputSource: make_msi_dir/directory

  tmb_dir:
    type: Directory
    outputSource: make_tmb_dir/directory
