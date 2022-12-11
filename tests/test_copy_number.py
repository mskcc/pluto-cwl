#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the copy_number.cwl
"""
import os
import sys
import unittest

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto import (
    CWLFile, 
    PlutoTestCase,
    DATA_SETS, 
    PORTAL_CNA_FILE,
    OFile
)
sys.path.pop(0)

class TestCopyNumber(PlutoTestCase):
    cwl_file = CWLFile('copy_number.cwl')

    def test_run_copy_number_one_file(self):
        """
        Test that Facets geneLevel copy number analysis step runs as expected with a single input file
        """
        self.input = {
            "output_cna_filename": "data_CNA.txt",
            "targets_list" : {
                "class": "File",
                "path": DATA_SETS['Proj_08390_G']['targets_list'],
            },
            "hisens_cncfs": [
                {
                    "class": "File",
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.cncf.txt")
                }
            ],
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_cna_file': OFile(name='data_CNA.txt', size=87905, hash='7cc89d24556de93b9a409812317581e67e5df494', dir = output_dir),
            'output_cna_ascna_file': OFile(name='data_CNA.ascna.txt', size=6164, hash='452d5ddef12a44693d5a98a05f5d300801734cfe', dir = output_dir),
            'output_cna_scna_file': OFile(name='data_CNA.scna.txt', size=5463, hash='8bec923ab1d622b4cf38ae042ac2416725650aed', dir = output_dir),
        }
        self.assertCWLDictEqual(output_json, expected_output)

        # load the data_CNA.txt file
        expected_header_parts = ['Tumor_Sample_Barcode', 'Hugo_Symbol', 'tcn', 'lcn', 'cf', 'tcn.em', 'lcn.em', 'cf.em', 'chr', 'seg.start', 'seg.end', 'frac_elev_major_cn', 'Nprobes', 'WGD', 'purity', 'FACETS_CALL', 'ccs_filter', 'review', 'FACETS_CNA']
        self.assertHeaderEquals(expected_output['output_cna_file']['path'], expected_header_parts)

        expected_header_parts = ['Hugo_Symbol', 's_C_VJ7F47_P001_d']
        self.assertHeaderEquals(expected_output['output_cna_ascna_file']['path'], expected_header_parts)

        expected_header_parts = ['Hugo_Symbol', 's_C_VJ7F47_P001_d']
        self.assertHeaderEquals(expected_output['output_cna_scna_file']['path'], expected_header_parts)

    def test_run_copy_number_two_files(self):
        """
        Test that Facets geneLevel copy number analysis step runs as expected with two input files
        """
        self.input = {
            "portal_CNA_file": "data_CNA.txt",
            "targets_list" : {
                "class": "File",
                "path": DATA_SETS['Proj_08390_G']['targets_list'],
            },
            "hisens_cncfs": [
                {
                    "class": "File",
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.cncf.txt")
                },
                {
                    "class": "File",
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample9.rg.md.abra.printreads__Sample10.rg.md.abra.printreads_hisens.cncf.txt")
                }
            ],
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_cna_file': OFile(name='data_CNA.txt', size=143118, hash='6dfa53b8a0fad1156060476bcf445d959f0e6eb2', dir = output_dir),
            'output_cna_ascna_file': OFile(name='data_CNA.ascna.txt', size=8658, hash='3953e55b3db85b69209982211c53b9d8f049dc01', dir = output_dir),
            'output_cna_scna_file': OFile(name='data_CNA.scna.txt', size=6937, hash='9ddcee42cce0d49aec5745303be480b6c4ef0fe8', dir = output_dir),
        }
        self.maxDiff = None
        self.assertCWLDictEqual(output_json, expected_output)

if __name__ == "__main__":
    unittest.main()
