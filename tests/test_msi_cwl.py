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
from pluto.settings import DATA_SETS
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
              "d":    {"class": "File", "path": "/work/ci/resources/request_files/msisensor/b37_known_somatic_microsatellites.list" },
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
                'checksum': 'sha1$576fa40d02976be7527b394f482d16768874b174',
                'size': 61,
                'path': os.path.join(output_dir, "Sample24.Sample23")
                }
            }
        self.assertDictEqual(output_json, expected_output)





if __name__ == "__main__":
    unittest.main()