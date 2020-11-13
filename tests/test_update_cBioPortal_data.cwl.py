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
    from .tools import load_mutations, run_cwl, write_table
    from .settings import CWL_DIR

if __name__ == "__main__":
    from tools import load_mutations, run_cwl, write_table
    from settings import CWL_DIR

cwl_file = os.path.join(CWL_DIR, 'update_cBioPortal_data.cwl')

class TestUpdate_cBioPortal_dataCWL(unittest.TestCase):
    def test_add_af(self):
        """
        Test Update_cBioPortal_dataCW with tiny dataset
        """
        # maf_lines = [
        #     ['# comment 1'],
        #     ['# comment 2'],
        #     ['Hugo_Symbol', 't_depth', 't_alt_count','tcn','lcn','expected_alt_copies','ccf_expected_copies','ccf_expected_copies_lower','ccf_expected_copies_upper'],
        #     ['SUFU', '100', '75','4','0.127','0.615','0.375'],
        #     ['GOT1', '100', '1' ,'4','0.127','0.615','0.375'], # need to change the values
        #     ['SOX9', '100', '0' ,'4','0.127','0.615','0.375'],
        # ]

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
                print_command = True
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'output.maf'),
                    'basename': 'output.maf',
                    'class': 'File',
                    #'checksum': 'sha1$39de59ad5d736db692504012ce86d3395685112e',
                    #'size': 109,
                    'path': os.path.join(output_dir, 'output.maf')
                    }
                }

            print('#######')
            print(input_json)
            print('#######')

            self.assertDictEqual(output_json, expected_output)

            comments, mutations = load_mutations(output_json['output_file']['path'])

            expected_comments = ['# comment 1', '# comment 2']
            self.assertEqual(comments, expected_comments)

            expected_mutations = [
                {'Hugo_Symbol': 'SUFU', 't_depth': '100', 't_alt_count':'75', 't_af': '0.75','tcn':'2', 'lcn':'1' ,'expected_alt_copies':'4','ccf_expected_copies':"0.127",'ccf_expected_copies_lower':"0.615",'ccf_expected_copies_upper':"0.375",'ASCN.TOTAL_COPY_NUMBER': "2",'ASCN.MINOR_COPY_NUMBER': "1",'ASCN.EXPECTED_ALT_COPIES': "4","ASCN.CCF_EXPECTED_COPIES": "0.127","ASCN.CCF_EXPECTED_COPIES_LOWER": "0.615","ASCN.CCF_EXPECTED_COPIES_UPPER": "0.375","ASCN.ASCN_METHOD": "FACETS","ASCN.ASCN_INTEGER_COPY_NUMBER": 'NA' },
                {'Hugo_Symbol': 'GOT1', 't_depth': '100', 't_alt_count':'1','t_af': '0.75','tcn':'2', 'lcn':'1' ,'expected_alt_copies':'4','ccf_expected_copies':"0.127",'ccf_expected_copies_lower':"0.615",'ccf_expected_copies_upper':"0.375",'ASCN.TOTAL_COPY_NUMBER': "2",'ASCN.MINOR_COPY_NUMBER': "1",'ASCN.EXPECTED_ALT_COPIES': "4","ASCN.CCF_EXPECTED_COPIES": "0.127","ASCN.CCF_EXPECTED_COPIES_LOWER": "0.615","ASCN.CCF_EXPECTED_COPIES_UPPER": "0.375","ASCN.ASCN_METHOD": "FACETS","ASCN.ASCN_INTEGER_COPY_NUMBER": 'NA' },
                {'Hugo_Symbol': 'SOX9', 't_depth': '100', 't_alt_count':'0','t_af': '0.75','tcn':'2', 'lcn':'1' ,'expected_alt_copies':'4','ccf_expected_copies':"0.127",'ccf_expected_copies_lower':"0.615",'ccf_expected_copies_upper':"0.375",'ASCN.TOTAL_COPY_NUMBER': "2",'ASCN.MINOR_COPY_NUMBER': "1",'ASCN.EXPECTED_ALT_COPIES': "4","ASCN.CCF_EXPECTED_COPIES": "0.127","ASCN.CCF_EXPECTED_COPIES_LOWER": "0.615","ASCN.CCF_EXPECTED_COPIES_UPPER": "0.375","ASCN.ASCN_METHOD": "FACETS","ASCN.ASCN_INTEGER_COPY_NUMBER": 'NA' }
                ]
            self.assertEqual(mutations, expected_mutations)


if __name__ == "__main__":
    unittest.main()
