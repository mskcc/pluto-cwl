#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
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
  mutation_maf_files:
    type: File[]
    doc: "analysis_mutations_filename (ANALYSIS_MUTATIONS_FILENAME) cbio_mutation_data_filename (CBIO_MUTATION_DATA_FILENAME): (MAF_DIR)/*.muts.maf"
  facets_hisens_seg_files:
    type: File[]
    doc: "cbio_segment_data_filename (CBIO_SEGMENT_DATA_FILENAME; <project_id>_data_cna_hg19.seg) analysis_segment_cna_filename (ANALYSIS_SEGMENT_CNA_FILE; <project_id>.seg.cna.txt): (FACETS_DIR)/*_hisens.seg"
  facets_hisens_cncf_files:
    type: File[]
    doc: "cbio_cna_data_filename (CBIO_CNA_DATA_FILENAME; data_CNA.txt) analysis_gene_cna_filename (ANALYSIS_GENE_CNA_FILENAME; <project_id>.gene.cna.txt): (FACETS_DIR)/*_hisens.cncf.txt"
  mutation_svs_maf_files:
    type: File[]
    doc: "analysis_sv_filename (ANALYSIS_SV_FILE; <project_id>.svs.maf): (MAF_DIR)/*.svs.pass.vep.maf"
  targets_list:
    type: File
  argos_version_string:
    type: string
    doc: "the version label of Roslin / Argos used to run the project analysis (ARGOS_VERSION_STRING)"
  is_impact:
    default: "True"
    type: string
    doc: "whether or not the project is an IMPACT project; should be the value 'True' if so, otherwise any other value means 'False' (IS_IMPACT)"
  helix_filter_version:
    type: string
    doc: "the version label of this helix filter repo (HELIX_FILTER_VERSION; git describe --all --long)"

steps:
  # <project_id>.gene.cna.txt (analysis_gene_cna_filename)
  generate_cna_data:
    run: copy_number.cwl
    in:
      output_cna_filename: analysis_gene_cna_filename
      output_cna_ascna_filename:
        valueFrom: ${ return inputs.output_cna_filename.replace(/\.[^/.]+$/, "") + '.ascna.txt'; }
      output_cna_scna_filename:
        valueFrom: ${ return inputs.output_cna_filename.replace(/\.[^/.]+$/, "") + '.scna.txt'; }
      targets_list: targets_list
      hisens_cncfs: facets_hisens_cncf_files
    out:
      [ output_cna_file ]

  # <project_id>.muts.maf (analysis_mutations_filename)
  # filter each maf file
  muts_maf_filter:
    run: maf_filter.cwl
    scatter: maf_file
    in:
      maf_file: mutation_maf_files
      argos_version_string: argos_version_string
      is_impact: is_impact
      analysis_mutations_filename: analysis_mutations_filename # <project_id>.muts.maf
    out: [ analysis_mutations_file ]
    # concat all the maf files into a single table
  concat_analysis_muts_maf:
    run: concat_with_comments.cwl
    in:
      input_files: muts_maf_filter/analysis_mutations_file
      comment_value: helix_filter_version
      output_filename: analysis_mutations_filename # <project_id>.muts.maf
    out:
      [ output_file ]

  # <project_id>.seg.cna.txt (analysis_segment_cna_filename)
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
  # rename the output file
  rename_analysis_hisens_segs:
    run: cp.cwl
    in:
      input_file: concat_hisens_segs/output_file
      output_filename: analysis_segment_cna_filename # <project_id>.seg.cna.txt
    out:
      [output_file]

  # <project_id>.svs.maf (analysis_sv_filename)
  # (MAF_DIR)/*.svs.pass.vep.maf (mutation_svs_maf_files)
  generate_analysis_svs_maf:
    run: concat_with_comments.cwl
    in:
      input_files: mutation_svs_maf_files
      comment_value: helix_filter_version
    out:
      [output_file]
  rename_analysis_svs_maf:
    run: cp.cwl
    in:
      input_file: generate_analysis_svs_maf/output_file
      output_filename: analysis_sv_filename # <project_id>.svs.maf
    out:
      [output_file]

  # create the 'analysis' directory and put some files in it
  make_analysis_dir:
    run: put_in_dir.cwl
    in:
      gene_cna_file: generate_cna_data/output_cna_file # <project_id>.gene.cna.txt
      muts_maf_file: concat_analysis_muts_maf/output_file # <project_id>.muts.maf
      hisens_segs: rename_analysis_hisens_segs/output_file # <project_id>.seg.cna.txt
      svs_maf_file: rename_analysis_svs_maf/output_file # <project_id>.svs.maf
      output_directory_name:
        valueFrom: ${ return "analysis"; }
      files:
        valueFrom: ${ return [
          inputs.gene_cna_file,
          inputs.muts_maf_file,
          inputs.hisens_segs,
          inputs.svs_maf_file
          ]}
    out: [ directory ]

outputs:
  analysis_dir:
    type: Directory
    outputSource: make_analysis_dir/directory
