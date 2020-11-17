#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the add_af.cwl
"""
import os
import unittest
from tempfile import TemporaryDirectory

# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import load_mutations, run_cwl, write_table, dicts2lines
    from .settings import CWL_DIR

if __name__ == "__main__":
    from tools import load_mutations, run_cwl, write_table, dicts2lines
    from settings import CWL_DIR

cwl_file = os.path.join(CWL_DIR, 'update_cBioPortal_data.cwl')

class TestUpdate_cBioPortal_dataCWL(unittest.TestCase):

    def setUp(self):
        self.maf_row1 = {
        "Hugo_Symbol" : "FGF3",
        "Entrez_Gene_Id" : "2248",
        "Chromosome" : "11",
        "Start_Position" : "69625447",
        "End_Position": "69625448",
        "Tumor_Sample_Barcode": "Sample1-T",
        "Matched_Norm_Sample_Barcode": "Sample1-N",
        "portal_val": "foo" # dummy value that would only be in portal data_mutations_extended.txt output maf file
        }
        self.maf_row2 = {
        "Hugo_Symbol" : "PNISR",
        "Entrez_Gene_Id" : "25957",
        "Chromosome" : "6",
        "Start_Position" : "99865784",
        "End_Position": "99865785",
        "Tumor_Sample_Barcode": "Sample1-T",
        "Matched_Norm_Sample_Barcode": "Sample1-N",
        "portal_val": "foo" # dummy value that would only be in portal data_mutations_extended.txt output maf file
        }
        self.maf_row3 = { # extra row with no match in facets
        "Hugo_Symbol" : "PNISR",
        "Entrez_Gene_Id" : "25957",
        "Chromosome" : "6",
        "Start_Position" : "99865788",
        "End_Position": "99865789",
        "Tumor_Sample_Barcode": "Sample1-T",
        "Matched_Norm_Sample_Barcode": "Sample1-N",
        "portal_val": "foo" # dummy value that would only be in portal data_mutations_extended.txt output maf file
        }
        self.facets_row1 = {
        "Hugo_Symbol" : "FGF3",
        "Entrez_Gene_Id" : "2248",
        "Chromosome" : "11",
        "Start_Position" : "69625447",
        "End_Position": "69625448",
        "Tumor_Sample_Barcode": "Sample1-T",
        "Matched_Norm_Sample_Barcode": "Sample1-N",
        "ASCN.TOTAL_COPY_NUMBER": "1" # dummy value that would only be in facets Tumor1.Normal1_hisens.ccf.maf output maf file
        }
        self.facets_row2 = {
        "Hugo_Symbol" : "PNISR",
        "Entrez_Gene_Id" : "25957",
        "Chromosome" : "6",
        "Start_Position" : "99865784",
        "End_Position": "99865785",
        "Tumor_Sample_Barcode": "Sample1-T",
        "Matched_Norm_Sample_Barcode": "Sample1-N",
        "ASCN.TOTAL_COPY_NUMBER": "2" # dummy value that would only be in facets Tumor1.Normal1_hisens.ccf.maf output maf file
        }
        self.facets_row3 = { # extra row with no match in maf
        "Hugo_Symbol" : "PNISR2",
        "Entrez_Gene_Id" : "25957",
        "Chromosome" : "6",
        "Start_Position" : "99865784",
        "End_Position": "99865785",
        "Tumor_Sample_Barcode": "Sample1-T",
        "Matched_Norm_Sample_Barcode": "Sample1-N",
        "ASCN.TOTAL_COPY_NUMBER": "2" # dummy value that would only be in facets Tumor1.Normal1_hisens.ccf.maf output maf file
        }

        self.demo_comments = [
        ['# comment 1'],
        ['# comment 2']
        ]


    def test_add_af(self):
        """
        Test Update_cBioPortal_dataCW with tiny dataset
        """

        self.maxDiff = None
        # make sets of lines to write to tables
        maf_rows = [ self.maf_row1, self.maf_row2, self.maf_row3 ]
        maf_lines = dicts2lines(dict_list = maf_rows, comment_list = self.demo_comments)

        facets_rows = [ self.facets_row1, self.facets_row2, self.facets_row3 ]
        facets_lines = dicts2lines(dict_list = facets_rows, comment_list = self.demo_comments)

        with TemporaryDirectory() as tmpdir:
            input_maf = write_table(tmpdir = tmpdir, filename = 'input.maf', lines = maf_lines)
            input_facets_file = write_table(tmpdir, filename = "facets.maf", lines = facets_lines)
            input_json = {
                "subcommand": "merge_mafs",
                "input_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "facets_maf":{
                      "class": "File",
                      "path": input_facets_file
                    },
                "output_filename":  'output.maf',
                }
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file,
                print_command = False,
                )

            expected_output = {
                'failed_txt': None,
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'output.maf'),
                    'basename': 'output.maf',
                    'class': 'File',
                    'checksum': 'sha1$73273a9ac68aa9a1b6e8d8747f268ea009bf7f0d',
                    'size': 335,
                    'path': os.path.join(output_dir, 'output.maf')
                    },
                'stderr_txt': {
                    'location': 'file://' + os.path.join(output_dir, 'output.maf_stderr.txt'),
                    'basename': 'output.maf_stderr.txt',
                    'class': 'File',
                    'checksum': 'sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709',
                    'size': 0,
                    'path': os.path.join(output_dir, 'output.maf_stderr.txt')
                    },
                'stdout_txt':{
                    'location': 'file://' + os.path.join(output_dir, 'output.maf_stdout.txt'),
                    'basename': 'output.maf_stdout.txt',
                    'class': 'File',
                    'checksum': 'sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709',
                    'size': 0,
                    'path': os.path.join(output_dir, 'output.maf_stdout.txt')}

                }

            self.assertDictEqual(output_json, expected_output)

            comments, mutations = load_mutations(output_json['output_file']['path'])

            expected_comments = ['# comment 1', '# comment 2']
            self.assertEqual(comments, expected_comments)

            expected_mutations = [
                {
                "Hugo_Symbol" : "FGF3",
                "Entrez_Gene_Id" : "2248",
                "Chromosome" : "11",
                "Start_Position" : "69625447",
                "End_Position": "69625448",
                "Tumor_Sample_Barcode": "Sample1-T",
                "Matched_Norm_Sample_Barcode": "Sample1-N",
                "portal_val": "foo",
                "ASCN.CLONAL": "1"
                },
                {
                "Hugo_Symbol" : "PNISR",
                "Entrez_Gene_Id" : "25957",
                "Chromosome" : "6",
                "Start_Position" : "99865784",
                "End_Position": "99865785",
                "Tumor_Sample_Barcode": "Sample1-T",
                "Matched_Norm_Sample_Barcode": "Sample1-N",
                "portal_val": "foo",
                "ASCN.CLONAL": "2"
                },
                {
                "Hugo_Symbol" : "PNISR",
                "Entrez_Gene_Id" : "25957",
                "Chromosome" : "6",
                "Start_Position" : "99865788",
                "End_Position": "99865789",
                "Tumor_Sample_Barcode": "Sample1-T",
                "Matched_Norm_Sample_Barcode": "Sample1-N",
                "portal_val": "foo",
                "ASCN.CLONAL": "."
                }
            ]
            self.assertEqual(mutations, expected_mutations)

if __name__ == "__main__":
    unittest.main()
