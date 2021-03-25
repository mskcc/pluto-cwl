#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the MSI analysis cwl
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
from pluto.settings import DATA_SETS, MICROSATELLITES_LIST
sys.path.pop(0)

class TestTMBWorkflow(PlutoTestCase):
    cwl_file = CWLFile('msi.cwl')
    def test_tmb_workflow(self):
        """
        Test case for the MSI analysis workflow
        """
        self.maxDiff = None

        output_file = os.path.join(self.tmpdir, "output.txt")

        self.input = {
              "d":    {"class": "File", "path": MICROSATELLITES_LIST },
              "n":    {"class": "File", "path": os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam") },
              "t":    {"class": "File", "path": os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam") },
              "o": "Sample24.Sample23"
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir, "Sample24.Sample23"),
                'basename': "Sample24.Sample23",
                'class': 'File',
                'checksum': 'sha1$8d5856202b859fc0da427bf9069eb870f2adcd55',
                'size': 62,
                'path': os.path.join(output_dir, "Sample24.Sample23")
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        lines = self.read_table(output_file)

        expected_lines = [
            ['Total_Number_of_Sites', 'Number_of_Somatic_Sites', '%'],
            ['628', '138', '21.97']
            ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
