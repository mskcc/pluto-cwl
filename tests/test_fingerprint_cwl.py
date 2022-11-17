#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the fingerprint.cwl
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
from pluto.serializer import OFile
from pluto.settings import CONPAIR_MARKERS_BED, CONPAIR_MARKERS_TXT
sys.path.pop(0)


class TestFingerprint(PlutoTestCase):
    cwl_file = CWLFile('fingerprint.cwl')

    def test_fingerprint_1(self):
        """
        """
        self.input = {
            "conpair_markers_txt": {
                "path": CONPAIR_MARKERS_TXT,
                "class": "File"
                },
            "conpair_markers_bed": {
                "path": CONPAIR_MARKERS_BED,
                "class": "File"
                },
            "ref_fasta": {
                "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA'],
                "class": "File"
                },
            "tumor_bam": {
                "path": os.path.join(self.DATA_SETS["demo"]["BAM_DIR"], "Sample23.bam"), #"/work/ci/dmp_finderprint_matching/dummy_bam/dummy.bam",
                "class": "File"
                },
            "dmp_dir": { ##"/work/ci/dmp/likelihoods"
                "class": "Directory",
                "path": "/work/ci/dmp_finderprint_matching/dummy_pickle"
                }
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name = 'Sample23.concordance.tsv', dir = output_dir)
            }

        # file contenst are inconsistent so strip some keys from the output dict
        strip_related_keys = [
            ('basename', 'Sample23.concordance.tsv', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)
        self.assertNumMutations(os.path.join(output_dir, 'Sample23.concordance.tsv'), 2)
        self.assertHeaderEquals(os.path.join(output_dir, 'Sample23.concordance.tsv'),
            ['concordance', 'num_markers_used', 'num_total_markers', "tumor", "normal", "tumor_pileup", "normal_pileup"])

        lines = self.read_table(os.path.join(output_dir, 'Sample23.concordance.tsv'))
        expected_lines = [
            ['concordance', 'num_markers_used', 'num_total_markers', 'tumor', 'normal'],
            ['0.4258517034068136', '998', '1024', 'Sample23', 'dummy'],
            ['0.4470468431771894', '982', '1024', 'Sample23', 'dummy']
            ]
        self.assertEqual([ l[0:5] for l in lines ], expected_lines)




if __name__ == "__main__":
    unittest.main()
