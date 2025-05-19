#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: ExpressionTool
doc: separate the samples list into samples with and without .maf files

requirements:
  - class: InlineJavascriptRequirement
  - $import: types.yml
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 3

inputs:
  samples:
    type: "types.yml#FilloutMafOptionalSample[]"

outputs:
  samples_with_maf: "types.yml#FilloutSample[]"
  samples_without_maf: "types.yml#FilloutNoMafsample[]"

expression: |
  ${
    var samples_with_maf = [];
    var samples_without_maf = [];

    for ( var i in inputs.samples ){
      var sample = inputs.samples[i];

      if (sample["maf_file"] == null ) {
        samples_without_maf.push(sample);
      } else {
        samples_with_maf.push(sample);
      };

    };

    var res = {
        "samples_with_maf": samples_with_maf,
        "samples_without_maf": samples_without_maf
      };
    return res;
  }

