#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
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
    default: true
    type: boolean
    doc: "whether or not the project is an IMPACT project; should be the value 'True' if so, otherwise any other value means 'False'"
  # TODO: this shouild actually be type: string[]
  extra_pi_groups:
    type: [ "null", string]
    default: null
    doc: "a list of other groups to be associated with the project in cBioPortal"
  cbio_segment_data_filename:
    type: string
    doc: "<project_id>_data_cna_hg19.seg)"
  cbio_meta_cna_segments_filename:
    type: string
    doc: "<project_id>_meta_cna_hg19_seg.txt)"
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
  cbio_meta_sv_filename:
    type: string
    default: meta_SV.txt
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
  cbio_sv_data_filename:
    type: string
    default: data_SV.txt
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
  mutation_maf_files:
    type: File[]
    doc: "*.muts.maf"
  facets_hisens_seg_files:
    type: File[]
    doc: "*_hisens.seg"
  facets_hisens_cncf_files:
    type: File[]
    doc: "*_hisens.cncf.txt"
  mutation_svs_txt_files:
    type: File[]
    doc: "*.svs.pass.vep.portal.txt"
  msi_files:
    type: File[]
    doc: "msi.tsv files"
  tmb_files:
    type: File[]
    doc: "*.tmb.tsv files"
  facets_suite_txt_files:
    type:
      - "null"
      - File[]
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
  update_extra_sample_ids:
    doc: if no extra sample ids were passed in, convert the null value to an empty list to make downstream processes easier
    in:
      sample_ids: extra_sample_ids
    out: [ extra_sample_ids ]
    run:
      class: ExpressionTool
      inputs:
        sample_ids:
          type:
          - "null"
          - string[]
      outputs:
        extra_sample_ids:
          type: string[]
      expression: |
        ${
          var sample_ids = [];
          if (inputs.sample_ids === null) {
            return {'extra_sample_ids': sample_ids};
          } else {
            return {'extra_sample_ids': inputs.sample_ids};
          }
        }

  # meta_clinical_sample.txt (cbio_clinical_sample_meta_filename; meta_clinical_sample_file)
  generate_meta_clinical_sample:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "meta_sample" }
      cancer_study_id: cancer_study_identifier
      sample_data_filename:  cbio_clinical_sample_data_filename # data_clinical_sample.txt
      output_filename: cbio_clinical_sample_meta_filename
    out:
      [output_file]

  # data_clinical_patient.txt (cbio_clinical_patient_data_filename; data_clinical_patient_file)
  generate_data_clinical_patient:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "patient" }
      data_clinical_file: data_clinical_file
      output_filename: cbio_clinical_patient_data_filename
    out:
      [output_file]

  # data_clinical_sample.txt (cbio_clinical_sample_data_filename)
  generate_data_clinical_sample:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "sample" }
      data_clinical_file: data_clinical_file
      sample_summary_file: sample_summary_file
      output_filename: cbio_clinical_sample_data_filename
      project_pi: project_pi
      request_pi: request_pi
      facets_txt_files: facets_suite_txt_files
    out:
      [output_file]

  # meta_study.txt (cbio_meta_study_filename; cbio_meta_study_file)
  generate_cbio_meta_study:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "study" }
      output_filename: cbio_meta_study_filename
      cancer_study_id: cancer_study_identifier
      name: project_name
      short_name: project_short_name
      type_of_cancer: cancer_type
      description: project_description
      extra_groups: extra_pi_groups
    out:
      [output_file]

  # meta_clinical_patient.txt (cbio_clinical_patient_meta_filename)
  generate_cbio_clinical_patient_meta:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "meta_patient" }
      output_filename: cbio_clinical_patient_meta_filename
      cancer_study_id: cancer_study_identifier
      patient_data_filename: cbio_clinical_patient_data_filename # data_clinical_patient.txt
    out:
      [output_file]

  # meta_CNA.txt (cbio_meta_cna_filename)
  generate_cbio_meta_cna:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "meta_cna" }
      output_filename: cbio_meta_cna_filename
      cancer_study_id: cancer_study_identifier
      cna_data_filename: cbio_cna_data_filename # data_CNA.txt
    out:
      [output_file]

  # meta_fusions.txt (cbio_meta_fusions_filename)
  generate_cbio_meta_fusions:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "meta_fusion" }
      output_filename: cbio_meta_fusions_filename
      cancer_study_id: cancer_study_identifier
      fusion_data_filename: cbio_fusion_data_filename # data_fusions.txt
    out:
      [output_file]

  # meta_SV.txt (cbio_meta_sv_filename)
  generate_cbio_meta_sv:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "meta_sv" }
      output_filename: cbio_meta_sv_filename
      cancer_study_id: cancer_study_identifier
      sv_data_filename: cbio_sv_data_filename # data_SV.txt
    out:
      [output_file]

  # meta_mutations_extended.txt (cbio_meta_mutations_filename)
  generate_meta_mutations_extended:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "meta_mutations" }
      output_filename: cbio_meta_mutations_filename
      cancer_study_id: cancer_study_identifier
      mutations_data_filename: cbio_mutation_data_filename # data_mutations_extended.txt
    out:
      [output_file]

  # <project_id>_meta_cna_hg19_seg.txt (cbio_meta_cna_segments_filename)
  generate_meta_cna_segments:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "meta_segments" }
      output_filename: cbio_meta_cna_segments_filename
      cancer_study_id: cancer_study_identifier
      segmented_data_filename: cbio_segment_data_filename # <project_id>_data_cna_hg19.seg
    out:
      [output_file]

  # cases_all.txt (cbio_cases_all_filename)
  generate_cbio_cases_all:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "cases_all" }
      output_filename: cbio_cases_all_filename
      cancer_study_id: cancer_study_identifier
      data_clinical_file: data_clinical_file
    out:
      [output_file]
  update_cases_all:
    run: updateCaseList.cwl
    in:
      sample_ids: update_extra_sample_ids/extra_sample_ids
      case_list: generate_cbio_cases_all/output_file
      output_filename: cbio_cases_all_filename
    out: [output_file]

  # cases_cnaseq.txt
  generate_cases_cnaseq:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "cases_cnaseq" }
      output_filename: cbio_cases_cnaseq_filename
      cancer_study_id: cancer_study_identifier
      data_clinical_file: data_clinical_file
    out:
      [output_file]
  update_cases_cnaseq:
    run: updateCaseList.cwl
    in:
      sample_ids: update_extra_sample_ids/extra_sample_ids
      case_list: generate_cases_cnaseq/output_file
      output_filename: cbio_cases_cnaseq_filename
    out: [output_file]

  # cases_cna.txt
  generate_cases_cna:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "cases_cna" }
      output_filename: cbio_cases_cna_filename
      cancer_study_id: cancer_study_identifier
      data_clinical_file: data_clinical_file
    out:
      [output_file]
  update_cases_cna:
    run: updateCaseList.cwl
    in:
      sample_ids: update_extra_sample_ids/extra_sample_ids
      case_list: generate_cases_cna/output_file
      output_filename: cbio_cases_cna_filename
    out: [output_file]

  # cases_sequenced.txt (cbio_cases_sequenced_filename)
  generate_cases_sequenced:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "cases_sequenced" }
      output_filename: cbio_cases_sequenced_filename
      cancer_study_id: cancer_study_identifier
      data_clinical_file: data_clinical_file
    out:
      [output_file]
  update_cases_sequenced:
    run: updateCaseList.cwl
    in:
      sample_ids: update_extra_sample_ids/extra_sample_ids
      case_list: generate_cases_sequenced/output_file
      output_filename: cbio_cases_sequenced_filename
    out: [output_file]


  # data_CNA.txt (cbio_cna_data_filename)
  # data_CNA.ascna.txt (cbio_cna_ascna_data_filename)
  # data_CNA.scna.txt, (cbio_cna_scna_data_filename)
  # (FACETS_DIR)/*_hisens.cncf.txt (facets_hisens_cncf_files)
  # targets_list
  generate_cna_data:
    run: copy_number.cwl
    in:
      output_cna_filename: cbio_cna_data_filename
      output_cna_ascna_filename: cbio_cna_ascna_data_filename
      output_cna_scna_filename: cbio_cna_scna_data_filename
      targets_list: targets_list
      hisens_cncfs: facets_hisens_cncf_files
    out:
      [ output_cna_file, output_cna_ascna_file, output_cna_scna_file ]
  # replace the 'ILLOGICAL' values in the data_CNA.scna.txt file
  # and output it as 'data_CNA.txt' instead
  replace_illogical_values:
    run: replace.cwl
    in:
      input_file: generate_cna_data/output_cna_scna_file # data_CNA.scna.txt
      output_filename: cbio_cna_data_filename # data_CNA.txt
    out:
      [ output_file ]
  # need to clean the header columns on some of the data_CNA.scna.txt and data_CNA.txt files
  clean_cna_headers:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "clean_cna" }
      input_file:  replace_illogical_values/output_file
      output_filename: cbio_cna_data_filename
    out:
      [output_file]
  clean_ascna_headers:
    run: generate_cBioPortal_file.cwl
    in:
      subcommand:
        valueFrom: ${ return "clean_cna" }
      input_file:  generate_cna_data/output_cna_ascna_file
      output_filename: cbio_cna_ascna_data_filename  # data_CNA.ascna.txt
    out:
      [output_file]
  # if there was extra CNA file, merge it in
  merge_cna:
    run: full-outer-join.cwl
    in:
      table1: clean_cna_headers/output_file
      table2: extra_cna_files
      join_key:
        valueFrom: ${ return "Hugo_Symbol" }
      output_filename:
        valueFrom: ${ return "data_CNA.txt" } # data_CNA_merged.txt
    out:
      [ output_file ]



  # data_mutations_extended.txt (cbio_mutation_data_filename)
  # filter each maf file
  muts_maf_filter:
    run: maf_filter.cwl
    scatter: maf_file
    in:
      maf_file: mutation_maf_files
      argos_version_string: argos_version_string
      is_impact: is_impact
      cbio_mutation_data_filename: cbio_mutation_data_filename # data_mutations_extended.txt
    out: [ cbio_mutation_data_file ]
    # concat all the maf files into a single table
  concat_cbio_muts_maf:
    run: concat-tables.cwl
    in:
      input_files: muts_maf_filter/cbio_mutation_data_file
      output_filename: cbio_mutation_data_filename # data_mutations_extended.txt
      comments:
        valueFrom: ${ return true; }
    out:
      [output_file]

  # <project_id>_data_cna_hg19.seg (cbio_segment_data_filename)
  # need to reduce the number of significant figures in the hisens_segs files
  reduce_sig_figs_hisens_segs:
    run: reduce_sig_figs.cwl
    scatter: input_file
    in:
      input_file: facets_hisens_seg_files
    out:
      [output_file]
  # concatenate all of the hisens_segs files
  concat_hisens_segs:
    run: concat.cwl
    in:
      input_files: reduce_sig_figs_hisens_segs/output_file
    out:
      [output_file]
  # rename the hisens_segs concatenated table to something that cBioPortal recognizes
  rename_cbio_hisens_segs:
    run: cp.cwl
    in:
      input_file: concat_hisens_segs/output_file
      output_filename: cbio_segment_data_filename # <project_id>_data_cna_hg19.seg
    out:
      [output_file]

  # data_fusions.txt (cbio_fusion_data_filename)
  # (mutation_svs_txt_files; (MAF_DIR)/*.svs.pass.vep.portal.txt)
  # concatenate all the mutation svs files
  generate_cbio_fusions_data:
    run: concat.cwl
    in:
      input_files: mutation_svs_txt_files
    out:
      [output_file]
  filter_cbio_fusions:
    run: fusion_filter.cwl
    in:
      fusions_file: generate_cbio_fusions_data/output_file
      output_filename: cbio_fusion_data_filename # data_fusions.txt
      known_fusions_file: known_fusions_file
    out:
      [output_file]
  convert_fusion_to_sv:
    run: fusion_to_sv.cwl
    in:
      fusion_file: filter_cbio_fusions/output_file
      output_filename: cbio_sv_data_filename # data_SV.txt
    out:
      [output_file]
  # create a case_list directory
  make_case_list_dir:
    run: put_in_dir.cwl
    in:
      # cases_all: generate_cbio_cases_all/output_file
      # cases_cnaseq: generate_cases_cnaseq/output_file
      # cases_cna: generate_cases_cna/output_file
      # cases_sequenced: generate_cases_sequenced/output_file
      cases_all: update_cases_all/output_file
      cases_cnaseq: update_cases_cnaseq/output_file
      cases_cna: update_cases_cna/output_file
      cases_sequenced: update_cases_sequenced/output_file
      output_directory_name:
        valueFrom: ${ return "case_lists"; }
      files:
        valueFrom: ${return [
          inputs.cases_all,
          inputs.cases_cnaseq,
          inputs.cases_cna,
          inputs.cases_sequenced
          ]}
    out: [ directory ]


#  HANDLING FOR MSI AND TMB FILES
  concat_tmb_tables:
    run: concat-tables.cwl
    in:
      input_files: tmb_files
      output_filename:
        valueFrom: ${ return "tmb.tsv"; }
      comments:
        valueFrom: ${ return true; }
    out: [ output_file ]


  concat_msi_tables:
    run: concat-tables.cwl
    in:
      input_files: msi_files
      output_filename:
        valueFrom: ${ return "msi.tsv"; }
      comments:
        valueFrom: ${ return true; }
    out:
      [ output_file ]


  # combine the MSI results with the data clinical file
  merge_msi_data_clinical:
    run: merge-tables.cwl
    in:
      table1: generate_data_clinical_sample/output_file
      table2: concat_msi_tables/output_file
      key1:
        valueFrom: ${ return "SAMPLE_ID"; } # sample column header from data clinical file
      key2:
        valueFrom: ${ return "SAMPLE_ID"; } # sample column header from MSI file
      output_filename:
        valueFrom: ${ return "data_clinical_sample.txt"; } # TODO: should this be passed in?
      cBioPortal:
        valueFrom: ${ return true; }
    out: [ output_file ]

  # # combine the TMB, MSI results with the data clinical file
  merge_tmb_data_clinical:
    run: merge-tables.cwl
    in:
      table1: merge_msi_data_clinical/output_file
      table2: concat_tmb_tables/output_file
      key1:
        valueFrom: ${ return "SAMPLE_ID"; } # sample column header from data clinical file
      key2:
        valueFrom: ${ return "SampleID"; } # sample column header from TMB file # TODO: This should be changed to SAMPLE_ID ?
      output_filename:
        valueFrom: ${ return "data_clinical_sample.txt"; } # TODO: should this be passed in?
      cBioPortal:
        valueFrom: ${ return true; }
    out: [ output_file ]




  # TODO: move this to workflow_with_facets.cwl
  compile_report:
    run: report.cwl
    in:
      mutation_file: concat_cbio_muts_maf/output_file
      samples_file: merge_tmb_data_clinical/output_file
      patients_file: generate_data_clinical_patient/output_file
    out: [ output_file ]




outputs:
  portal_meta_clinical_sample_file:
    type: File
    outputSource: generate_meta_clinical_sample/output_file # meta_clinical_sample.txt
  portal_data_clinical_patient_file:
    type: File
    outputSource: generate_data_clinical_patient/output_file # data_clinical_patient.txt
  portal_data_clinical_sample_file:
    type: File
    outputSource: merge_tmb_data_clinical/output_file # data_clinical_sample.txt
  portal_meta_study_file:
    type: File
    outputSource: generate_cbio_meta_study/output_file # meta_study.txt
  portal_clinical_patient_meta_file:
    type: File
    outputSource: generate_cbio_clinical_patient_meta/output_file # meta_clinical_patient.txt
  portal_meta_cna_file:
    type: File
    outputSource: generate_cbio_meta_cna/output_file # meta_CNA.txt
  portal_meta_fusions_file:
    type: File
    outputSource: generate_cbio_meta_fusions/output_file # meta_fusions.txt
  portal_meta_sv_file:
    type: File
    outputSource: generate_cbio_meta_sv/output_file # meta_SV.txt
  portal_meta_mutations_extended_file:
    type: File
    outputSource: generate_meta_mutations_extended/output_file # meta_mutations_extended.txt
  portal_meta_cna_segments_file:
    type: File
    outputSource: generate_meta_cna_segments/output_file  # <project_id>_meta_cna_hg19_seg.txt
  portal_cna_data_file:
    type: File
    # outputSource: clean_cna_headers/output_file # data_CNA.txt
    outputSource: merge_cna/output_file # data_CNA.txt
  portal_cna_ascna_file:
    type: File
    outputSource: clean_ascna_headers/output_file # data_CNA.ascna.txt
  portal_muts_file:
    type: File
    outputSource: concat_cbio_muts_maf/output_file # data_mutations_extended.txt
  portal_hisens_segs:
    type: File
    outputSource: rename_cbio_hisens_segs/output_file # # <project_id>_data_cna_hg19.seg
  portal_fusions_data_file:
    type: File
    outputSource: filter_cbio_fusions/output_file # data_fusions.txt
  portal_sv_data_file:
    type: File
    outputSource: convert_fusion_to_sv/output_file # data_SV.txt
  portal_case_list_dir:
    type: Directory
    outputSource: make_case_list_dir/directory
  portal_report:
    type: File
    outputSource: compile_report/output_file
