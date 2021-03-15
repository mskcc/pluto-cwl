#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
doc: "
Special case of concat-tables.py to use with concatenating maf files for
use with https://github.com/zengzheng123/GetBaseCountsMultiSample

Note that these columns are required by GetBaseCountsMultiSample:
( as per https://github.com/zengzheng123/GetBaseCountsMultiSample/blob/master/GetBaseCountsMultiSample.cpp#L771-L801 )

Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, t_ref_count, t_alt_count, n_ref_count, n_alt_count, Variant_Classification

However, only these columns appear to influence the output:

Hugo_Symbol Chromosome Start_Position End_Position Variant_Classification Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2

usually in a column order matching this;
grep -v '#' Proj_08390_G.muts.maf | cut -f1,5,6,7,9,11,12,13,16,17,41,42,44,45 > Proj_08390_G.cols.muts.maf

also, if present these columns also appear in output but seem to be optional

Center	NCBI_Build

see also;
https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
"

baseCommand: [
  "concat-tables.py",
  '--dir',
  '--na-str', '.', # empty string here crashes GetBaseCountsMultiSample
  # keep only these cols in the output;
  '--keep-cols', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Variant_Classification', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Center', 'NCBI_Build',
  # replace the value of these cols with the na_str to save space in the output file;
  '--na-cols', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 't_ref_count', 't_alt_count', 'n_ref_count', 'n_alt_count',
  '--comments', '--no-carriage-returns'
  ]

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: mskcc/helix_filters_01:21.02.2
  InitialWorkDirRequirement:
    listing:
      - entryname: inputs_dir
        writable: true
        # put all the input files in a dir and rename them so the filenames do not collide
        entry: |-
            ${
              for (var i = 0; i < inputs.input_files.length; i++) {
                inputs.input_files[i].basename = inputs.input_files[i].basename + "." + i;
              }
              return {class: 'Directory', listing: inputs.input_files};
            }

arguments:
  - valueFrom: $(inputs.output_filename)
    position: 1
    prefix: -o
  - valueFrom: ${ return "inputs_dir" }
    position: 2

inputs:
  output_filename:
    type: string
    default: "output.maf"
  input_files:
    type: File[]

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
