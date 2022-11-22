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
            "dmp_dir": {
                "class": "Directory",
                "path": self.DATA_SETS["Conpair_1"]["LIKELIHOODS"]
                },
            "additional_normal_bams": [
                {"class": "File", "path": os.path.join(self.DATA_SETS["Conpair_1"]["BAM_DIR"], "Sample1.UnitTest01.bam")},
                {"class": "File", "path": os.path.join(self.DATA_SETS["demo"]["BAM_DIR"], "Sample23.bam")}, # test against itself
            ]
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name = 'Sample23.concordance.tsv', dir = output_dir, hash = "dd689be28d551b19033db5e120286d1e1e8d83f2", size = 692)
            }

        lines = self.read_table(os.path.join(output_dir, 'Sample23.concordance.tsv'))
        expected_lines = [
            ['concordance', 'num_markers_used', 'num_total_markers', 'tumor', 'normal', 'tumor_filename', 'normal_filename'],
            ['0.35070140280561124', '998', '1024', 'Sample23', 'Sample5', 'Sample23.pileup', 'Sample5.UnitTest01.pickle'],
            ['0.34468937875751504', '998', '1024', 'Sample23', 'Sample1-foo', 'Sample23.pileup', 'Sample1.UnitTest01.pickle'],
            ['0.3486973947895792', '998', '1024', 'Sample23', 'Sample2-bar', 'Sample23.pileup', 'Sample2.UnitTest01.pickle'],
            ['0.3408856848609681', '971', '1024', 'Sample23', 'Sample3-baz', 'Sample23.pileup', 'Sample3.UnitTest01.pickle'],
            ['0.35015447991761073', '971', '1024', 'Sample23', 'Sample4', 'Sample23.pileup', 'Sample4.UnitTest01.pickle'],
            ['1.0', '999', '1024', 'Sample23', 'Sample23', 'Sample23.pileup', 'Sample23.pileup'],
            ['0.34468937875751504', '998', '1024', 'Sample23', 'Sample1', 'Sample23.pileup', 'Sample1.UnitTest01.pileup']
            ]

        self.assertEqual(lines, expected_lines)

        self.assertCWLDictEqual(output_json, expected_output)
        self.assertNumMutations(os.path.join(output_dir, 'Sample23.concordance.tsv'), 7)
        self.assertHeaderEquals(os.path.join(output_dir, 'Sample23.concordance.tsv'),
            ['concordance', 'num_markers_used', 'num_total_markers', "tumor", "normal", "tumor_filename", "normal_filename"])






if __name__ == "__main__":
    unittest.main()
