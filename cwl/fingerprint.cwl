#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
doc: "
Workflow to run generate fingerprint search on entire dmp normals and optionally on project normals
"

requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}

# NOTE: the following section is needed because the input dmp_dir has 10's of thousands of files and that makes the CWL blow up unless you turn off this setting
hints:
  cwltool:LoadListingRequirement:
    loadListing: no_listing
$namespaces:
  cwltool: "http://commonwl.org/cwltool#"

inputs:
  dmp_dir: # TODO: rename this to something more generic since we might not always be using DMP files
    doc: directory of reference samples to run concordance against
    type: Directory
    default: # TODO: don't put a default here, instead make this an optional "none"
      class: Directory
      path: /work/ci/dmp_finderprint_matching/dummy_pileup
  conpair_markers_bed:
    type: File
  conpair_markers_txt:
    type: File
  tumor_bam: # GATK .pileup or likelihoods .pickle
    # TODO: make this a list of samples, each with a bam
    type: File
    secondaryFiles: ["^.bai"]
  additional_normal_bams:
    doc: list of individual extra sample files to include with concordance
    type: File[]
    secondaryFiles: ["^.bai"]
    default: [{class: File, path: /work/ci/dmp_finderprint_matching/dummy_bam/dummy.bam}]
    # TODO: remove this default
  ref_fasta:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
      - .fai
      - ^.dict

outputs:
  output_file:
    type: File
    outputSource: run_conpair_concordance/output_file

steps:
  run_pileup_tumor:
    doc: create a pileup file for the input tumor samples
    run: conpair-pileup.cwl
    in:
      bam: tumor_bam
      ref: ref_fasta
      gatk:
        valueFrom: ${ return '/usr/bin/gatk.jar'; }
      markers_bed: conpair_markers_bed
      java_xmx:
        valueFrom: ${ return ["24g"]; }
      outfile:
        valueFrom: ${ return inputs.bam.basename.replace(".bam", ".pileup"); }
    out: [out_file]



  run_pileup_additional_normals:
    doc: create pileup files for the extra normal files provided
    run: conpair-pileup.cwl
    scatter: bam
    in:
        bam: additional_normal_bams
        ref: ref_fasta
        gatk:
          valueFrom: ${ return '/usr/bin/gatk.jar'; }
        markers_bed: conpair_markers_bed
        java_xmx:
          valueFrom: ${ return ["24g"]; }
        outfile:
          valueFrom: ${ return inputs.bam.basename.replace(".bam", ".pileup"); }
    out: [out_file]


  put_in_dir:
    doc: put the pileups for the extra normals into a dir for easier handling
    run: put_in_dir.cwl
    in:
      output_directory_name:
        valueFrom: ${ return "additional_normals"; }
      files: run_pileup_additional_normals/out_file
    out: [ directory ]


  run_conpair_concordance:
    doc: run Conpair on the set of tumor and normal files
    in:
      dmp_dir: dmp_dir
      markers: conpair_markers_txt
      tumor_file: run_pileup_tumor/out_file
      additional_normal_pickles: put_in_dir/directory
    out: [output_file]
    run:
      class: CommandLineTool
      baseCommand: ['bash', 'run.sh']
      requirements:
        ResourceRequirement:
          # ramMin: 16000
          coresMin: 8 # NOTE: make sure this matches what is used in the command
        DockerRequirement:
          dockerPull: mskcc/conpair:1.0.1
        InitialWorkDirRequirement:
          listing:
            - entryname: run.sh
              entry: |-
                set -eu
                # need to make a list for all the input files
                # NOTE: this needs to be done here because the files get staged in a container with tmp paths
                path="${ return inputs.dmp_dir.path; }"
                find "\$path"/ -type f > normal_dmp_pickles_file_list.txt

                path="${ return inputs.additional_normal_pickles.path; }"
                find "\$path"/ -type f >> normal_dmp_pickles_file_list.txt

                tumor_pickle="${ return inputs.tumor_file.path; }"
                markers="${ return inputs.markers.path; }"
                run.py \\
                  concordance \\
                  \$tumor_pickle \\
                  --normals-list normal_dmp_pickles_file_list.txt \\
                  --markers \$markers \\
                  --threads 8 \\
                  --output-file "${ return inputs.tumor_file.nameroot + '.concordance.tsv' }"
      inputs:
        dmp_dir: Directory
        additional_normal_pickles: Directory
        markers: File
        tumor_file: File

      outputs:
        output_file:
          type: File
          outputBinding:
            glob: ${ return inputs.tumor_file.nameroot + '.concordance.tsv' }



