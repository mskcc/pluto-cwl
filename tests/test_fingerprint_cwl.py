#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the paste-col.cwl
"""
import os
import sys
import unittest
# import json
# from tempfile import TemporaryDirectory

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

        self.preserve = True

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
            "runparams": {
                "tmp_dir": "/scratch",
                "gatk_jar_path": "/usr/bin/gatk.jar"
                },
            "tumor_bam": {
                "path": "/work/ci/dmp_finderprint_matching/dummy_bam/dummy.bam",
                "class": "File"
                },
            "dmp_dir": {
                "class": "Directory",
                "path": "/work/ci/dmp_finderprint_matching/dummy_pickle"
                }
            }



##"/work/ci/dmp/likelihoods"

        output_json, output_dir = self.run_cwl()
        
        expected_output = {
            'output_file': OFile(name = 'dummy.concordance.tsv', hash = "sha1$6b1a5ef498be9c15c842033ba7b5cc970397b5f5", size = 487, dir = output_dir)
            }


        self.assertCWLDictEqual(expected_output, output_json)




if __name__ == "__main__":
    unittest.main()
