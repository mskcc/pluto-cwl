#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
doc: index a Fillout sample .bam file
baseCommand: [ "bash", "run.index_sample_bam.sh" ]
requirements:
  - $import: types.yml
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 4 # make sure this matches what is in the command below!
  - class: DockerRequirement
    dockerPull: mskcc/helix_filters_01:samtools-1.9
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.sample['bam_file']) # NOTE: why are we doing this?? # NOTE: I think this is causing all input bam files to get copied when using Toil??
      - entryname: run.index_sample_bam.sh
        entry: |-
          set -eu
          # sample.bam
          input_bam="$(inputs.sample['bam_file'].basename)"
          # sample.bam.bai
          default_bai="\${input_bam}.bai"
          # sample.bai
          extra_bai="\${input_bam%.*}.bai"
          samtools index -@ 4 "\${input_bam}"
          cp "\${default_bai}" "\${extra_bai}"

inputs:
  sample: "types.yml#FilloutIndexSample"
outputs:
  sample:
    type: "types.yml#FilloutIndexedSample"
    outputBinding:
      outputEval: ${
        var ret = inputs.sample;
        ret['bam_file']['secondaryFiles'] = [{"class":"File", "path":runtime.outdir + "/" + inputs.sample["bam_file"].nameroot + ".bai"}];
        return ret;
        }
