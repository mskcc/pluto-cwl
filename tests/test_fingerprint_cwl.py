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
sys.path.pop(0)


class TestFingerprint(PlutoTestCase):
    cwl_file = CWLFile('fingerprint.cwl')

    def test_fingerprint_1(self):
        """
        """
        self.input = {
            "conpair_markers_txt": {
                "path": "/juno/work/ci/concordance-workflow/markers/IMPACT468/FP_tiling_genotypes_for_Conpair.txt",
                "class": "File"
                },
            "conpair_markers_bed": {
                "path": "/juno/work/ci/concordance-workflow/markers/IMPACT468/FP_tiling_genotypes_for_Conpair.bed",
                "class": "File"
                },
            "ref_fasta": {
                "path": "/juno/work/ci/resources/genomes/GRCh37/fasta/b37.fasta",
                "class": "File"
                },
            "tumor_bam": {
                "path": "/work/ci/dmp_finderprint_matching/dummy_bam/dummy.bam",
                "class": "File"
                },
            "dmp_dir": { ##"/work/ci/dmp/likelihoods"
                "class": "Directory",
                "path": "/work/ci/dmp_finderprint_matching/dummy_pickle"
                }
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name = 'dummy.concordance.tsv', dir = output_dir)
            }

        # file contenst are inconsistent so strip some keys from the output dict
        strip_related_keys = [
            ('basename', 'dummy.concordance.tsv', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)
        self.assertNumMutations(os.path.join(output_dir, 'dummy.concordance.tsv'), 2)
        self.assertHeaderEquals(os.path.join(output_dir, 'dummy.concordance.tsv'),
            ['concordance', 'num_markers_used', 'num_total_markers', "tumor", "normal", "tumor_pileup", "normal_pileup"])




if __name__ == "__main__":
    unittest.main()
