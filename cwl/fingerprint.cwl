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

hints:
  cwltool:LoadListingRequirement:
    loadListing: no_listing

inputs:
  dmp_dir:
    type: Directory
    default:
      class: Directory
      path: /work/ci/dmp_finderprint_matching/dummy_pileup
  conpair_markers_bed:  # conpair_markers_txt
    type: File
  conpair_markers_txt:  # conpair_markers_txt
    type: File
  tumor_bam: # GATK .pileup or likelihoods .pickle
    type: File
    secondaryFiles: ["^.bai"]
  additional_normal_bams:
    type: File[]
    secondaryFiles: ["^.bai"]
    default: [{class: File, path: /work/ci/dmp_finderprint_matching/dummy_bam/dummy.bam}]
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
  run-pileup-tumor:
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
######
  run-pileup-additional_normals:
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

  put-in-dir:
    run: put_in_dir.cwl
    in:
      output_directory_name:
        valueFrom: ${ return "additional_normals"; }
      files: run-pileup-additional_normals/out_file
    out: [ directory ]
#######
  run_conpair_concordance:
    in:
      dmp_dir: dmp_dir
      markers: conpair_markers_txt
      tumor_file: run-pileup-tumor/out_file
      additional_normal_pickles: put-in-dir/directory

    out: [output_file]

    run:
      class: CommandLineTool
      baseCommand: ['bash', 'run.sh']
      requirements:
        DockerRequirement:
          dockerPull: mskcc/conpair:1.0.1
        InitialWorkDirRequirement:
          listing:
            # NOTE: might need dos2unix for some that give errors ERROR: Your MAF uses CR line breaks, which we can't support. Please use LF or CRLF.
            # NOTE: might also need sanity check that maf has >1 line
            - entryname: run.sh
              entry: |-
                set -eu
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
                  --threads 10 \\
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



$namespaces:
  cwltool: "http://commonwl.org/cwltool#"
