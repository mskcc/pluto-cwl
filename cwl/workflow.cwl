#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: "
CWL workflow for generating Roslin / Argos post pipeline analysis files and cBioPortal data and metadata files

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
normal_bam_files
tumor_bam_files
assay_coverage

The following filenames are required:

analysis_mutations_filename
analysis_gene_cna_filename
analysis_sv_filename
analysis_segment_cna_filename
cbio_segment_data_filename
cbio_meta_cna_segments_filename
microsatellites_file

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
    ├── data_fusions.txt
    ├── data_mutations_extended.txt
    ├── meta_clinical_patient.txt
    ├── meta_clinical_sample.txt
    ├── meta_CNA.txt
    ├── meta_fusions.txt
    ├── meta_mutations_extended.txt
    ├── meta_study.txt
    ├── <project_id>_data_cna_hg19.seg
    └── <project_id>_meta_cna_hg19_seg.txt
"

requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  project_id:
    type: string
    doc: "unique identifier for the project (PROJ_ID)"
  project_pi:
    type: string
    doc: "principle investigator for the project (PROJ_PI)"
  request_pi:
    type: string
    doc: "principle investigator who requested the project (REQUEST_PI)"
  project_short_name:
    type: string
    doc: "a short name for the project in cBioPortal (PROJ_SHORT_NAME)"
  project_name:
    type: string
    doc: "a formal name for the project (PROJ_NAME)"
  project_description:
    type: string
    doc: "a description of the project (PROJ_DESC)"
  cancer_type:
    type: string
    doc: "the type of cancer used in the project (CANCER_TYPE)"
  cancer_study_identifier:
    type: string
    doc: "a study identifier for the project to use in cBioPortal (CANCER_STUDY_IDENTIFIER)"
  argos_version_string:
    type: string
    doc: "the version label of Roslin / Argos used to run the project analysis (ARGOS_VERSION_STRING)"
  helix_filter_version:
    type: string
    doc: "the version label of this helix filter repo (HELIX_FILTER_VERSION; git describe --all --long)"
  is_impact:
    default: True
    type: boolean
    doc: "whether or not the project is an IMPACT project; should be the value 'True' if so, otherwise any other value means 'False' (IS_IMPACT)"
  # TODO: this shouild actually be type: string[]
  extra_pi_groups:
    type: [ "null", string]
    default: null
    doc: "a list of other groups to be associated with the project in cBioPortal (EXTRA_PI_GROUPS)"
  analysis_segment_cna_filename:
    type: string
    doc: "(ANALYSIS_SEGMENT_CNA_FILE; <project_id>.seg.cna.txt)"
  analysis_sv_filename:
    type: string
    doc: "(ANALYSIS_SV_FILE; <project_id>.svs.maf)"
  analysis_gene_cna_filename:
    type: string
    doc: "(ANALYSIS_GENE_CNA_FILENAME; <project_id>.gene.cna.txt)"
  analysis_mutations_filename:
    type: string
    doc: "(ANALYSIS_MUTATIONS_FILENAME; <project_id>.muts.maf)"
  analysis_mutations_share_filename:
    type: string
    doc: "<project_id>.muts.share.maf)"
  cbio_segment_data_filename:
    type: string
    doc: "(CBIO_SEGMENT_DATA_FILENAME; <project_id>_data_cna_hg19.seg)"
  cbio_meta_cna_segments_filename:
    type: string
    doc: "(cbio_meta_cna_segments_filename; <project_id>_meta_cna_hg19_seg.txt)"
  cbio_cases_sequenced_filename:
    type: string
    doc: "(CBIO_CASES_SEQUENCED_FILE)"
    default: cases_sequenced.txt
  cbio_cases_cna_filename:
    type: string
    default: cases_cna.txt
    doc: "(CBIO_CASES_CNA_FILE)"
  cbio_cases_cnaseq_filename:
    type: string
    default: cases_cnaseq.txt
    doc: "(CBIO_CASES_CNASEQ_FILE)"
  cbio_cases_all_filename:
    type: string
    default: cases_all.txt
    doc: "(CBIO_CASES_ALL_FILE)"
  cbio_meta_mutations_filename:
    type: string
    default: meta_mutations_extended.txt
    doc: "(CBIO_META_MUTATIONS_FILE)"
  cbio_meta_fusions_filename:
    type: string
    default: meta_fusions.txt
    doc: "(CBIO_META_FUSIONS_FILE)"
  cbio_meta_cna_filename:
    type: string
    default: meta_CNA.txt
    doc: "(CBIO_META_CNA_FILE)"
  cbio_meta_study_filename:
    type: string
    default: meta_study.txt
    doc: "(CBIO_META_STUDY_FILE)"
  cbio_clinical_patient_meta_filename:
    type: string
    default: meta_clinical_patient.txt
    doc: "(CBIO_CLINCAL_PATIENT_META_FILE)"
  cbio_clinical_sample_meta_filename:
    type: string
    default: meta_clinical_sample.txt
    doc: "(CBIO_CLINICAL_SAMPLE_META_FILE)"
  cbio_clinical_sample_data_filename:
    type: string
    default: data_clinical_sample.txt
    doc: "(CBIO_CLINICAL_SAMPLE_DATA_FILENAME)"
  cbio_clinical_patient_data_filename:
    type: string
    default: data_clinical_patient.txt
    doc: "(CBIO_CLINCIAL_PATIENT_DATA_FILENAME)"
  cbio_fusion_data_filename:
    type: string
    default: data_fusions.txt
    doc: "(CBIO_FUSION_DATA_FILENAME)"
  cbio_mutation_data_filename:
    type: string
    default: data_mutations_extended.txt
    doc: "(CBIO_MUTATION_DATA_FILENAME)"
  cbio_cna_data_filename:
    type: string
    default: data_CNA.txt
    doc: "(CBIO_CNA_DATA_FILENAME)"
  cbio_cna_ascna_data_filename:
    type: string
    default: data_CNA.ascna.txt
    doc: "(CBIO_CNA_ASCNA_DATA_FILE)"
  cbio_cna_scna_data_filename:
    type: string
    default: data_CNA.scna.txt
    doc: "(CBIO_CNA_SCNA_DATA_FILE)"
  # mutation_maf_files:
  #   type: File[]
  #   doc: "analysis_mutations_filename (ANALYSIS_MUTATIONS_FILENAME) cbio_mutation_data_filename (CBIO_MUTATION_DATA_FILENAME): (MAF_DIR)/*.muts.maf"

  facets_hisens_seg_files:
    type:
      - "null"
      - File[]
      - string
    doc: "cbio_segment_data_filename (CBIO_SEGMENT_DATA_FILENAME; <project_id>_data_cna_hg19.seg) analysis_segment_cna_filename (ANALYSIS_SEGMENT_CNA_FILE; <project_id>.seg.cna.txt): (FACETS_DIR)/*_hisens.seg"
    default: [ {"path": "/work/ci/vurals/helix_manual_reruns/test_wo_facets_workflow/hisens.seg", "class": "File"} ]
      # valueFrom: ${ return {"path": "/work/ci/vurals/helix_manual_reruns/test_wo_facets_workflow/hisens.seg", "class": "File"}; }

  facets_hisens_cncf_files:
    type:
      - "null"
      - File[]
    doc: "cbio_cna_data_filename (CBIO_CNA_DATA_FILENAME; data_CNA.txt) analysis_gene_cna_filename (ANALYSIS_GENE_CNA_FILENAME; <project_id>.gene.cna.txt): (FACETS_DIR)/*_hisens.cncf.txt"
    default: [ {"path": "/work/ci/vurals/helix_manual_reruns/test_wo_facets_workflow/data_CNA.txt", "class": "File"} ]
      # valueFrom: ${ return [{"path": "/work/ci/vurals/helix_manual_reruns/test_wo_facets_workflow/data_CNA.txt", "class": "File"}]; }


  mutation_svs_txt_files:
    type: File[]
    doc: "cbio_fusion_data_filename (CBIO_FUSION_DATA_FILENAME; data_fusions.txt): (MAF_DIR)/*.svs.pass.vep.portal.txt"
  mutation_svs_maf_files:
    type: File[]
    doc: "analysis_sv_filename (ANALYSIS_SV_FILE; <project_id>.svs.maf): (MAF_DIR)/*.svs.pass.vep.maf"
  facets_suite_txt_files:
    type:
      - "null"
      - File[]
    doc: "Facets Suite .txt files for all samples in the request"
  # facets_aggregate_file:
  #   doc: "Facets Suite .txt file aggregated for all samples in the request"
  #   type:
  #     - "null"
  #     - File
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
  pairs:
    type:
      type: array
      items:
        - type: record
          fields:
            pair_maf: File
            snp_pileup: File
            pair_id: string
            tumor_id: string
            normal_id: string
  IMPACT_gene_list:
    type: File
    doc: "TSV file with gene labels and corresponding impact assays"

  assay_coverage:
    type: string
    doc: "genome_coverage value; amount of the genome in bp covered by the assay"

  microsatellites_file:
    type: File
    doc: "Microsatellites list file to use with MSI Sensor"

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




steps:
  run_analysis_workflow:
    run: analysis-workflow_wo_facets.cwl
    in:
      analysis_segment_cna_filename: analysis_segment_cna_filename
      analysis_sv_filename: analysis_sv_filename
      analysis_gene_cna_filename: analysis_gene_cna_filename
      analysis_mutations_filename: analysis_mutations_filename
      analysis_mutations_share_filename: analysis_mutations_share_filename
      pair: pairs
      mutation_maf_files:
        valueFrom: ${ return [inputs.pair[0].pair_maf] }


      facets_hisens_seg_files: facets_hisens_seg_files
          # valueFrom: ${ return [{"path": "/work/ci/vurals/helix_manual_reruns/test_wo_facets_workflow/hisens.seg", "class": "File"}]; }
      facets_hisens_cncf_files: facets_hisens_cncf_files
        # valueFrom: ${ return ["/work/ci/vurals/helix_manual_reruns/test_wo_facets_workflow/data_CNA.txt"]; }


      mutation_svs_maf_files: mutation_svs_maf_files
      targets_list: targets_list
      argos_version_string: argos_version_string
      is_impact: is_impact
      helix_filter_version: helix_filter_version
      IMPACT_gene_list: IMPACT_gene_list
    out:
      [ analysis_dir ]

  run_portal_workflow:
    run: portal-workflow_wo_facets.cwl
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

      pair: pairs
      mutation_maf_files:
        valueFrom: ${ return [inputs.pair[0].pair_maf] }


      facets_hisens_seg_files: facets_hisens_seg_files
        # valueFrom: ${ return [ {"path": "/work/ci/vurals/helix_manual_reruns/test_wo_facets_workflow/hisens.seg", "class": "File"}  ]; }
      facets_hisens_cncf_files: facets_hisens_cncf_files
        # valueFrom: ${ return ["/work/ci/vurals/helix_manual_reruns/test_wo_facets_workflow/data_CNA.txt"] }


      mutation_svs_txt_files: mutation_svs_txt_files
      targets_list: targets_list
      known_fusions_file: known_fusions_file
      data_clinical_file: data_clinical_file
      sample_summary_file: sample_summary_file
      facets_suite_txt_files: facets_suite_txt_files
    out:
      # [ portal_dir ]
      [
      portal_meta_clinical_sample_file, # meta_clinical_sample.txt
      portal_data_clinical_patient_file, # data_clinical_patient.txt
      portal_data_clinical_sample_file, # data_clinical_sample.txt
      portal_meta_study_file, # meta_study.txt
      portal_clinical_patient_meta_file, # meta_clinical_patient.txt
      # portal_meta_cna_file, # meta_CNA.txt
      portal_meta_fusions_file, # meta_fusions.txt
      portal_meta_mutations_extended_file, # meta_mutations_extended.txt
      # portal_meta_cna_segments_file, # <project_id>_meta_cna_hg19_seg.txt
      # portal_cna_data_file, # data_CNA.txt
      # portal_cna_ascna_file, # data_CNA.ascna.txt
      portal_muts_file, # data_mutations_extended.txt
      # portal_hisens_segs, # <project_id>_data_cna_hg19.seg
      portal_fusions_data_file, # data_fusions.txt
      portal_case_list_dir,
      # merged_cna_file, # data_CNA_merged.txt -> data_CNA.txt
      portal_report
      ]



  # run the TMB workflow
  run_tmb_workflow:
    run: tmb_workflow.cwl
    in:
      data_clinical_file: run_portal_workflow/portal_data_clinical_sample_file
      assay_coverage: assay_coverage
      pairs: pairs
    out:
      [ output_file ] # updated data_clinical_sample_file with the new TMB data

  # run the MSI workflow
  run_msi_workflow:
    run: msi_workflow.cwl
    in:
      data_clinical_file: run_portal_workflow/portal_data_clinical_sample_file # run_tmb_workflow/output_file # data_clinical_sample.txt
      microsatellites_file: microsatellites_file
      pairs: pairs
      normal_bam_files: normal_bam_files
      tumor_bam_files: tumor_bam_files
    out:
      [ output_file ] # updated data_clinical_file with MSI scores

  # combine the TMB, MSI results with the data clinical file
  merge_data_clinical:
    run: merge-tables.cwl
    in:
      table1: run_tmb_workflow/output_file
      table2: run_msi_workflow/output_file
      key1:
        valueFrom: ${ return "SAMPLE_ID"; } # sample column header from data clinical file
      key2:
        valueFrom: ${ return "SAMPLE_ID"; } # sample column header from MSI file
      output_filename:
        valueFrom: ${ return "data_clinical_sample.txt"; } # TODO: should this be passed in?
      cBioPortal:
        valueFrom: ${ return true; }
    out:
      [ output_file ]






  # create the "portal" directory in the output dir and put cBioPortal files in it
  make_portal_dir:
    run: put_in_dir.cwl
    in:
      portal_meta_clinical_sample_file: run_portal_workflow/portal_meta_clinical_sample_file # meta_clinical_sample.txt
      portal_data_clinical_patient_file: run_portal_workflow/portal_data_clinical_patient_file # data_clinical_patient.txt
      portal_data_clinical_sample_file: merge_data_clinical/output_file # data_clinical_sample.txt
      portal_meta_study_file: run_portal_workflow/portal_meta_study_file # meta_study.txt
      portal_clinical_patient_meta_file: run_portal_workflow/portal_clinical_patient_meta_file # meta_clinical_patient.txt
      # portal_meta_cna_file: run_portal_workflow/portal_meta_cna_file # meta_CNA.txt
      portal_meta_fusions_file: run_portal_workflow/portal_meta_fusions_file # meta_fusions.txt
      portal_meta_mutations_extended_file: run_portal_workflow/portal_meta_mutations_extended_file # meta_mutations_extended.txt
      # portal_meta_cna_segments_file: run_portal_workflow/portal_meta_cna_segments_file  # <project_id>_meta_cna_hg19_seg.txt
      # portal_cna_data_file: run_portal_workflow/merged_cna_file # data_CNA.txt
      # portal_cna_ascna_file: run_portal_workflow/portal_cna_ascna_file # data_CNA.ascna.txt
      portal_muts_file: run_portal_workflow/portal_muts_file # data_mutations_extended.txt
      # portal_hisens_segs: run_portal_workflow/portal_hisens_segs # # <project_id>_data_cna_hg19.seg
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
          inputs.portal_cna_data_file,
          inputs.portal_cna_ascna_file,
          inputs.portal_muts_file,
          inputs.portal_hisens_segs,
          inputs.portal_fusions_data_file,
          inputs.portal_case_list_dir,
          ]}
    out: [ directory ]

outputs:
  portal_dir:
    type: Directory
    outputSource: make_portal_dir/directory

  analysis_dir:
    type: Directory
    outputSource: run_analysis_workflow/analysis_dir
