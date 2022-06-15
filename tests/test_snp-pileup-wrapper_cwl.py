#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the snp-pileup-wrapper.cwl file
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
from pluto.settings import FACETS_SNPS_VCF, DATA_SETS
from pluto.serializer import OFile
sys.path.pop(0)

class TestSnpPileupCWL(PlutoTestCase):
    cwl_file = CWLFile('snp-pileup-wrapper.cwl')
    def test_snp_pileup1(self):
        """
        """
        self.input = {
            "snps_vcf": {
                "path": FACETS_SNPS_VCF,
                "class": "File"
            },
            # the smallest pair of bam files in the test dataset
            "normal_bam": {
                "path": os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam"),
                "class": "File"
            },
            "tumor_bam": {
                "path": os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam"),
                "class": "File"
            },
            "output_prefix": "Sample24.Sample23"
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name = "Sample24.Sample23.snp_pileup.gz", size = 34851004, hash = "755a8b64f45c819b4e2c481e64bf2fe36d1f5361", dir = output_dir)
            }
        self.maxDiff = None
        self.assertCWLDictEqual(output_json, expected_output)

if __name__ == "__main__":
    unittest.main()
