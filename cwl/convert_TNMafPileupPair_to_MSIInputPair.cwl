#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: ExpressionTool
requirements:
- $import: types.yml

inputs:
  pairs: "types.yml#TNMafPileupPair[]"

outputs:
  pairs: "types.yml#MSIInputPair[]"

expression: |
  ${
    var pairs = [];
    for ( var i in inputs.pairs ){
      var pair = {
        "tumor_id": inputs.pairs[i].tumor_id,
        "normal_id": inputs.pairs[i].normal_id,
        "pair_id": inputs.pairs[i].pair_id,
      };
      pairs.push(pair);
    };
    return {"pairs": pairs};
  }
