#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the maf_filter.cwl
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import CWLFile, PlutoTestCase
from pluto.settings import DATA_SETS, ARGOS_VERSION_STRING, IS_IMPACT
from pluto.serializer import OFile
sys.path.pop(0)

class TestMafFilter(PlutoTestCase):
    cwl_file = CWLFile('maf_filter.cwl')

    def test_filter_a_maf_file(self):
        """
        Test that a filtered maf file comes out as expected
        """
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        self.assertNumMutations(input_maf, 12514)

        self.input = {
            "maf_file": {
                  "class": "File",
                  "path": input_maf
                },
            "argos_version_string": ARGOS_VERSION_STRING,
            "is_impact": True,
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "cbio_mutation_data_filename": 'data_mutations_extended.txt'
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_mutations_file': OFile(name='Proj_08390_G.muts.maf', size=28079, hash='24421ab8d1a39a71f48eecbb0dd167d5d9f5c529', dir = output_dir),
            'cbio_mutation_data_file': OFile(name='data_mutations_extended.txt', size=4534, hash='6131494536ce956d741c820378e7e2ce1c714403', dir = output_dir),
            'rejected_file': OFile(name='rejected.muts.maf', size=18627626, hash='a06789623715703c5006db6876ecb58b8498f938', dir = output_dir),
            }
        self.assertCWLDictEqual(output_json, expected_output)

        self.assertNumMutations(os.path.join(output_dir, "Proj_08390_G.muts.maf"), 22)

        # validate output mutation file contents
        self.assertCompareMutFiles(
            os.path.join(output_dir, "Proj_08390_G.muts.maf"),
            os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Sample1", "analyst_file.txt"),
            muts_only = True,
            compare_len = True
        )

        self.assertCompareMutFiles(
            os.path.join(output_dir, 'data_mutations_extended.txt'),
            os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Sample1", "portal_file.txt"),
            muts_only = True,
            compare_len = True
        )

    def test_maf_filter_argos_3_2_0(self):
        """
        Test the maf filter script results when used with argos_version_string 3.2.0
        """
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        self.assertNumMutations(input_maf, 12514)

        self.input = {
            "maf_file": {
                  "class": "File",
                  "path": input_maf
                },
            "argos_version_string": "3.2.0",
            "is_impact": True,
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "cbio_mutation_data_filename": 'data_mutations_extended.txt'
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_mutations_file': OFile(name='Proj_08390_G.muts.maf', size=28081, hash='fd78842c9410e7e622dee270ec9c0e7628811f18', dir = output_dir),
            'cbio_mutation_data_file': OFile(name='data_mutations_extended.txt', size=4536, hash='47e716eabbfda3408b2d9a08b9bb432b2cb8fce8', dir = output_dir),
            'rejected_file': OFile(name='rejected.muts.maf', size=18627626, hash='a06789623715703c5006db6876ecb58b8498f938', dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)

        self.assertNumMutations(expected_output['analysis_mutations_file']['path'], 22)

        # validate output mutation file contents
        self.assertCompareMutFiles(
            expected_output['analysis_mutations_file']['path'],
            os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Sample1", "analyst_file.txt"),
            muts_only = True,
            compare_len = True
        )

        self.assertCompareMutFiles(
            expected_output['cbio_mutation_data_file']['path'],
            os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Sample1", "portal_file.txt"),
            muts_only = True,
            compare_len = True
        )


    def test_filter_maf_file_impact_false(self):
        """
        Test that a filtered maf file comes out as expected
        """
        self.maxDiff = None
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        self.assertNumMutations(input_maf, 12514)

        self.input = {
            "maf_file": {
                  "class": "File",
                  "path": input_maf
                },
            "argos_version_string": ARGOS_VERSION_STRING,
            "is_impact": False,
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "cbio_mutation_data_filename": 'data_mutations_extended.txt'
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_mutations_file': OFile(name='Proj_08390_G.muts.maf', size=24524, hash='9fb9d43c71e546750ddec6aea2313dda28547b3a', dir = output_dir),
            'cbio_mutation_data_file': OFile(name='data_mutations_extended.txt', size=3931, hash='15ca06249511c32c32e058c246a757ec8df11d83', dir = output_dir),
            'rejected_file': OFile(name='rejected.muts.maf', size=18790398, hash='e7441703699e82cef500d9557bfcbd3464ce8eab', dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)
        self.assertNumMutations(expected_output['analysis_mutations_file']['path'], 18)

    def test_large_maf_file(self):
        """
        Test that a giant maf file with tons of variants gets filtered as expected
        """
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Proj_08390_G", "Proj_08390_G.muts.maf")

        self.input = {
            "maf_file": {
                  "class": "File",
                  "path": input_maf
                },
            "argos_version_string": "2.x",
            "is_impact": True,
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "cbio_mutation_data_filename": 'data_mutations_extended.txt'
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_mutations_file': OFile(name='Proj_08390_G.muts.maf', size=2386906, hash='4ef341ab4280140f9be15e65a0258a4170ff651d', dir = output_dir),
            'cbio_mutation_data_file': OFile(name='data_mutations_extended.txt', size=278458, hash='af36cf815820fdf41f1401578138b5cbd551a217', dir = output_dir),
            'rejected_file': OFile(name='rejected.muts.maf', size=1047796463, hash='345953da2c7cb801fa08368260469cf7c153055f', dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)

        self.assertNumMutations(expected_output['analysis_mutations_file']['path'], 1662)
        self.assertNumMutations(expected_output['cbio_mutation_data_file']['path'], 1139)

        # validate output mutation file contents
        self.assertCompareMutFiles(
            expected_output['analysis_mutations_file']['path'],
            os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Proj_08390_G", "analyst_file.txt"),
            muts_only = True,
            compare_len = True
        )

        self.assertCompareMutFiles(
            expected_output['cbio_mutation_data_file']['path'],
            os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Proj_08390_G", "portal_file.txt"),
            muts_only = True,
            compare_len = True
        )

if __name__ == "__main__":
    unittest.main()
