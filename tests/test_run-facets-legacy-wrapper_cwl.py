#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the run-facets-wrapper.cwl file
"""
import os
import sys
import unittest



from pluto import (
    PlutoTestCase, 
    CWLFile,
    DATA_SETS,
    OFile
)


class TestRunFacetsLegacyWrapperCWL(PlutoTestCase):
    cwl_file = CWLFile('run-facets-legacy-wrapper.cwl')
    def test_run_facets_wrapper(self):
        """
        """
        self.input = {
            "snp_pileup": {
                "path": os.path.join(DATA_SETS['Proj_08390_G']['SNP_PILEUP_DIR'], "Sample24.Sample23.snp_pileup.gz"),
                "class": "File"
            },
            "sample_id": "Sample24"
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'arm_level_txt': OFile(name = "Sample24.arm_level.txt", hash = "df37c54ae4969257e436a7a7a595c42ef19ecbb5", size = 1824, dir = output_dir),
            'failed': False,
            'gene_level_txt': OFile(name = "Sample24.gene_level.txt", hash = "4e916a52458151007486bf536acfff539fdc2ecc", size = 148195, dir = output_dir),
            'hisens_Rdata': OFile(name = "Sample24_hisens.Rdata", hash = "6cc3fd1fad17111e32c7c88b259a523092539181", size = 214260, dir = output_dir),
            'hisens_seg': OFile(name = "Sample24_hisens.seg", hash = "652f9c6d0245af49bac6ca67a089af7d4e46801b", size = 1897, dir = output_dir),
            'purity_Rdata': OFile(name = "Sample24_purity.Rdata", hash = "ca2c46ceebbb960a02fb960134ba1711983e71c8", size = 213542, dir = output_dir),
            'purity_seg': OFile(name = "Sample24_purity.seg", hash = "591e6d8b432e1e910fe4fb4b1814508131f960c9", size = 1285, dir = output_dir),
            'qc_txt': OFile(name = "Sample24.qc.txt", hash = "d4a36726a5fcb7b268aae02d97ce4e382e42d9f6", size = 1339, dir = output_dir),
            'hisens_cncf_txt': OFile(name = "Sample24_hisens.cncf.txt", hash = "db9131a33889a1cac82e3bd6b3f0e5e182c65105", size = 5238, dir = output_dir),
            'purity_cncf_txt': OFile(name = "Sample24_purity.cncf.txt", hash = "b331530e1e46b5ba1bdcedeb67f2aa82da6ebc5f", size = 3630, dir = output_dir),
            'stderr_txt': OFile(name = "facets_legacy_stderr.txt", dir = output_dir),
            'stdout_txt': OFile(name = 'facets_legacy_stdout.txt', dir = output_dir),
        }
        self.maxDiff = None

        strip_related_keys = [
        ('basename', 'facets_legacy_stdout.txt', ['size', 'checksum']),
        ('basename', 'facets_legacy_stderr.txt', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        with open(os.path.join(output_dir, 'Sample24_hisens.cncf.txt')) as fin:
            self.assertEqual(len(fin.readlines()), 37)



