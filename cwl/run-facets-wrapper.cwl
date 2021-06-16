#!/usr/bin/env cwl-runner

# run-facets-wrapper.R \
# --counts-file $${snp_pileup} \
# --sample-id $${pair_id} \
# --purity-cval 100 \
# --cval 50 \
# --seed 1000 \
# --everything \
# --min-nhet 25 \
# --purity-min-nhet 25 \
# -D $(FACETS_OUTPUT) \
# --facets-lib-path /usr/local/lib/R/site-library

cwlVersion: v1.0
class: CommandLineTool
successCodes: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255]
stdout: facets_stdout.txt
stderr: facets_stderr.txt
baseCommand: ["bash", "run_facets_wrapper.sh"]

requirements:
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:facets-suite-2.0.6
  InitialWorkDirRequirement:
    listing:
      - entryname: run_facets_wrapper.sh
        entry: |-
          run-facets-wrapper.R --everything -D . --facets-lib-path /usr/local/lib/R/site-library $@ || touch failed.txt

inputs:
  snp_pileup:
    type: File
    inputBinding:
      position: 1
      prefix: '--counts-file'
  sample_id:
    type: string
    inputBinding:
      position: 2
      prefix: '--sample-id'
  purity_cval:
    type: ["null", string]
    default: "100"
    inputBinding:
      position: 3
      prefix: '--purity-cval'
  cval:
    type: ["null", string]
    default: "50"
    inputBinding:
      position: 4
      prefix: '--cval'
  seed:
    type: ["null", string]
    default: "1000"
    inputBinding:
      position: 5
      prefix: '--seed'
  min_nhet:
    type: ["null", string]
    default: "25"
    inputBinding:
      position: 6
      prefix: '--min-nhet'
  purity_min_nhet:
    type: ["null", string]
    default: "25"
    inputBinding:
      position: 7
      prefix: '--purity-min-nhet'

# s_C_ABCD_P001_d_s_C_ABCD_N001_d_purity.seg
# s_C_ABCD_P001_d_s_C_ABCD_N001_d_purity.png
# s_C_ABCD_P001_d_s_C_ABCD_N001_d_hisens.seg
# s_C_ABCD_P001_d_s_C_ABCD_N001_d_hisens.png
# s_C_ABCD_P001_d_s_C_ABCD_N001_d.qc.txt
# s_C_ABCD_P001_d_s_C_ABCD_N001_d.gene_level.txt
# s_C_ABCD_P001_d_s_C_ABCD_N001_d.arm_level.txt
# s_C_ABCD_P001_d_s_C_ABCD_N001_d.txt
# s_C_ABCD_P001_d_s_C_ABCD_N001_d_purity.rds
# s_C_ABCD_P001_d_s_C_ABCD_N001_d_hisens.rds
outputs:
  purity_seg:
    type: File?
    outputBinding:
      glob: $(inputs.sample_id)_purity.seg
  hisens_seg:
    type: File?
    outputBinding:
      glob: $(inputs.sample_id)_hisens.seg
  qc_txt:
    type: File?
    outputBinding:
      glob: $(inputs.sample_id).qc.txt
  gene_level_txt:
    type: File?
    outputBinding:
      glob: $(inputs.sample_id).gene_level.txt
  arm_level_txt:
    type: File?
    outputBinding:
      glob: $(inputs.sample_id).arm_level.txt
  output_txt:
    type: File?
    outputBinding:
      glob: $(inputs.sample_id).txt
  purity_rds:
    type: File?
    outputBinding:
      glob: $(inputs.sample_id)_purity.rds
  hisens_rds:
    type: File?
    outputBinding:
      glob: $(inputs.sample_id)_hisens.rds
  failed_txt:
    type: File?
    outputBinding:
      glob: failed.txt
  stdout_txt:
    type: stdout
  stderr_txt:
    type: stderr
