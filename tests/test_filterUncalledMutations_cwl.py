#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the filterUncalledMutations.cwl
"""
import os
import sys
import unittest

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto import (
    DATA_SETS,
    PlutoTestCase, 
    CWLFile,
    OFile, 
    ODir
)
sys.path.pop(0)


class TestFilterUncalledMutations(PlutoTestCase):
    cwl_file = CWLFile('filterUncalledMutations.cwl')

    def test_1(self):
        maf_file = os.path.join(DATA_SETS["Fillout01"]["OUTPUT_DIR"], "output.maf")

        comments, mutations = self.load_mutations(maf_file, strip = True)
        self.assertEqual(len(mutations), 475)

        self.input = {
        "input_file": {"class": "File", "path": maf_file}
        }

        output_json, output_dir = self.run_cwl()

        output_data_mutations_extended = os.path.join(output_dir,'data_mutations_extended.txt')
        output_data_mutations_uncalled = os.path.join(output_dir,'data_mutations_uncalled.txt')

        expected_output = {
            'called_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir, hash = 'e7430656d9fcbce36fa57eb92460db57742168ae', size = 347254),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir, hash = '58129786cc299011202eb078734b3ff513d54081', size = 287883),
        }

        self.maxDiff = None

        self.assertCWLDictEqual(output_json, expected_output)

        comments, mutations = self.load_mutations(output_data_mutations_extended, strip = True)
        self.assertEqual(len(mutations), 253)

        comments, mutations = self.load_mutations(output_data_mutations_uncalled, strip = True)
        self.assertEqual(len(mutations), 222)


if __name__ == "__main__":
    unittest.main()
