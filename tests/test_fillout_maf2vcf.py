#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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


class TestFilloutMaf2Vcf(PlutoTestCase):
    cwl_file = CWLFile('fillout_maf2vcf.cwl')

    def test_maf2vcf(self):
        """
        """
        sample1_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample1.FillOutUnitTest01.muts.maf')
        sample2_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample2.FillOutUnitTest01.muts.maf')
        sample3_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample3.FillOutUnitTest01.muts.maf')
        sample4_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample4.FillOutUnitTest01.muts.maf')
        sample5_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample5.FillOutUnitTest01.muts.maf')

        self.input = {
            "sample_id": "Sample1",
            "maf_file": {"class": "File", "path": sample1_maf},
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']}
        }
        output_json, output_dir = self.run_cwl()
        expected_output = {
            "output_file": OFile(name="Sample1.sorted.vcf.gz", 
                size = 2923, 
                hash = "73a623dd614467e9069ca5ba74b9da16cd881af0", 
                dir = output_dir, 
                secondaryFiles = [OFile(name = "Sample1.sorted.vcf.gz.tbi", 
                    hash = "f176cd2d3572eb0502637db147f82ca55a23c1fd", 
                    size = 3842, 
                    dir = output_dir)]),
            "output_vcf": OFile(name = "Sample1.sorted.vcf", hash = "26a3d2b3bda116805fdf1e4b820981b830497a82", size = 9956, dir = output_dir)
        }
        # fields inside the vcf are not static due to timestamps, etc..
        strip_related_keys = [
        ('basename', 'Sample1.sorted.vcf.gz', ['size', 'checksum']),
        ('basename', 'Sample1.sorted.vcf.gz.tbi', ['size', 'checksum']),
        ('basename', 'Sample1.sorted.vcf', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)
        self.assertNumMutationsHash(expected_output["output_file"].path, 65, "5a45a60bee13b211dd8c2b6082d7e83f")
        self.assertNumMutationsHash(expected_output["output_vcf"].path, 65, "5a45a60bee13b211dd8c2b6082d7e83f")
    




    # def test_convert_all(self):
    #     """
    #     Not a real test case... 
    #     """
    #     from pprint import pprint
    #     sample1_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample1.FillOutUnitTest01.muts.maf')
    #     sample2_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample2.FillOutUnitTest01.muts.maf')
    #     sample3_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample3.FillOutUnitTest01.muts.maf')
    #     sample4_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample4.FillOutUnitTest01.muts.maf')
    #     sample5_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample5.FillOutUnitTest01.muts.maf')
    #     samples = [
    #         ("Sample1", sample1_maf),
    #         ("Sample2", sample2_maf),
    #         ("Sample3", sample3_maf),
    #         ("Sample4", sample4_maf),
    #         ("Sample5", sample5_maf),
    #         ]
    #     for sample in samples:
    #         self.input = {
    #             "sample_id": sample[0],
    #             "maf_file": {"class": "File", "path": sample[1]},
    #             "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']}
    #         }
    #         output_json, output_dir = self.run_cwl()
    #         pprint(output_json)
            



