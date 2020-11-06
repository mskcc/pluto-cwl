#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the maf_col_filter.cwl
"""
import os
import unittest
from tempfile import TemporaryDirectory

# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import run_command, load_mutations, run_cwl, write_table
    from .settings import CWL_DIR, DATA_SETS

if __name__ == "__main__":
    from tools import run_command, load_mutations, run_cwl, write_table
    from settings import CWL_DIR, DATA_SETS

cwl_file = os.path.join(CWL_DIR, 'maf_col_filter.cwl')


# copy/pasted this from maf_col_filter.py for testing purposes
cols_to_keep = set([
"Hugo_Symbol",
"Entrez_Gene_Id",
"Center",
"NCBI_Build",
"Chromosome",
"Start_Position",
"End_Position",
"Strand",
"Variant_Classification",
"Variant_Type",
"Reference_Allele",
"Tumor_Seq_Allele1",
"Tumor_Seq_Allele2",
"dbSNP_RS",
"dbSNP_Val_Status",
"Tumor_Sample_Barcode",
"Matched_Norm_Sample_Barcode",
"Match_Norm_Seq_Allele1",
"Match_Norm_Seq_Allele2",
"Tumor_Validation_Allele1",
"Tumor_Validation_Allele2",
"Match_Norm_Validation_Allele1",
"Match_Norm_Validation_Allele2",
"Verification_Status",
"Validation_Status",
"Mutation_Status",
"Sequencing_Phase",
"Sequence_Source",
"Validation_Method",
"Score",
"BAM_File",
"Sequencer",
"Tumor_Sample_UUID",
"Matched_Norm_Sample_UUID",
"HGVSc",
"HGVSp",
"Amino_Acid_Change", # added by maf_filter script so not in raw maf
"Transcript_ID",
"Exon_Number",
"t_depth",
"t_ref_count",
"t_alt_count",
"n_depth",
"n_ref_count",
"n_alt_count"
])


class TestMafColFilter(unittest.TestCase):
    def test_filter_maf_file_cols(self):
        """
        Filter columns in a tiny demo maf file
        """
        maf_lines = [
            ['# comment 1'], # keep the comments
            ['# comment 2'],
            ['Hugo_Symbol', 'foo_value'], # foo_value column should be removed in output
            ['SUFU', '1'],
            ['GOT1', '2']
        ]
        # run the script in a temporary directory
        with TemporaryDirectory() as tmpdir:
            input_maf_file = write_table(tmpdir = tmpdir, filename = 'input.maf', lines = maf_lines)
            input_json = {
                "input_file": {
                      "class": "File",
                      "path": input_maf_file
                    },
                "output_filename": "output.maf"
            }

            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'output.maf'),
                    'basename': 'output.maf',
                    'class': 'File',
                    'checksum': 'sha1$e55f7bdaa146f37b48d6c920ed27184e394ef1e6',
                    'size': 46,
                    'path': os.path.join(output_dir, 'output.maf')
                    }
                }
            self.assertDictEqual(output_json, expected_output)

            # validate number of lines output
            with open(output_json['output_file']['path']) as fin:
                output_maf_lines = len(fin.readlines())
            self.assertEqual(output_maf_lines, 5)

            # validate file contents
            comments, mutations = load_mutations(output_json['output_file']['path'])

            expected_comments = ['# comment 1', '# comment 2']
            self.assertEqual(comments, expected_comments)

            expected_mutations = [{'Hugo_Symbol': 'SUFU'}, {'Hugo_Symbol': 'GOT1'}]
            self.assertEqual(mutations, expected_mutations)

    def test_filter_maf_file_cols_full(self):
        """
        Test col filter on a full sized dataset
        """
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")

        with TemporaryDirectory() as tmpdir:
            input_json = {
                "input_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "output_filename": "output.maf"
            }

            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'output.maf'),
                    'basename': 'output.maf',
                    'class': 'File',
                    'checksum': 'sha1$a2f5b9f1533fd443b41561ca718ffca62ab45f36',
                    'size': 2710681,
                    'path': os.path.join(output_dir, 'output.maf')
                    }
                }
            self.assertDictEqual(output_json, expected_output)

            # validate number of lines output
            with open(output_json['output_file']['path']) as fin:
                output_maf_lines = len(fin.readlines())
            self.assertEqual(output_maf_lines, 12518)

            # validate file contents
            comments, mutations = load_mutations(output_json['output_file']['path'])

            self.assertEqual(len(mutations), 12514)

            for key in mutations[0].keys():
                self.assertTrue(key in cols_to_keep)

            # make sure there are fewer than or equal to the number of columns in new output as there are entries to keep 
            self.assertTrue( len(mutations[0].keys()) <= len(cols_to_keep) )


if __name__ == "__main__":
    unittest.main()
