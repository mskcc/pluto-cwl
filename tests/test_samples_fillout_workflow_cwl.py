#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the samples_fillout_workflow cwl
"""
import os
import sys
import unittest
from collections import OrderedDict

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile, TableReader, md5_obj
from pluto.settings import ENABLE_LARGE_TESTS, DATA_SETS
from pluto.serializer import OFile, ODir
sys.path.pop(0)

from fixtures_fillout import rows

class TestSamplesFillout(PlutoTestCase):
    cwl_file = CWLFile('samples_fillout_workflow.cwl')

    def test_Nick_testcase(self):
        """
        Test case using Nick's custom made maf and bam files for fillout testing

        This test cases uses the germline filter to exclude some mutations in the output

        Takes about 10min to run
        """
        self.maxDiff = None
        self.runner_args['use_cache'] = False # do not use cache because it breaks for some reason
        self.runner_args['debug'] = True
        self.runner_args['js_console'] = True

        sample1_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample1.FillOutUnitTest01.muts.maf')
        sample2_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample2.FillOutUnitTest01.muts.maf')
        sample3_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample3.FillOutUnitTest01.muts.maf')
        sample4_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample4.FillOutUnitTest01.muts.maf')
        sample5_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample5.FillOutUnitTest01.muts.maf')

        sample1_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample1.UnitTest01.bam')
        sample2_bam =os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample2.UnitTest01.bam')
        sample3_bam =os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample3.UnitTest01.bam')
        sample4_bam =os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample4.UnitTest01.bam')
        sample5_bam =os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample5.UnitTest01.bam')

        self.input = {
            "samples": [
                {
                    "sample_id": "Sample1",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample1_maf },
                    "bam_file": { "class": "File", "path": sample1_bam }
                },
                {
                    "sample_id": "Sample2",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample2_maf },
                    "bam_file": { "class": "File", "path": sample2_bam }
                },
                {
                    "sample_id": "Sample3",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "clinical",
                    "maf_file": { "class": "File", "path": sample3_maf },
                    "bam_file": { "class": "File", "path": sample3_bam }
                },
                {
                    "sample_id": "Sample4",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "clinical",
                    "maf_file": { "class": "File", "path": sample4_maf },
                    "bam_file": { "class": "File", "path": sample4_bam }
                },
                {
                    "sample_id": "Sample5",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample5_maf },
                    "bam_file": { "class": "File", "path": sample5_bam }
                },
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']}
        }
        output_json, output_dir = self.run_cwl()
        output_path = os.path.join(output_dir,'output.maf')
        filtered_output_path = os.path.join(output_dir,'output.filtered.maf')
        portal_output_path = os.path.join(output_dir,'data_mutations_extended.txt')
        uncalled_output_path = os.path.join(output_dir,'data_mutations_uncalled.txt')

        expected_output = {
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        }

        # file contenst are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)
        # all_effects field is variable and changes bytes and checksum
        # need to check number of variant outputs instead
        self.assertNumMutationsHash(output_path, 475, 'd041bc641d85761b60c6b7ef8606bab2')
        self.assertNumMutationsHash(filtered_output_path, 230, 'c9cde01507d1b2470057c5d120eaab68')
        self.assertNumMutationsHash(portal_output_path, 163, '8dd6f3af030a2eca3b5fa0698896361a')
        self.assertNumMutationsHash(uncalled_output_path, 67, 'a474b61268d2a4c25fd27cc2ccbbce96')
        self.assertEqualNumMutations([portal_output_path, uncalled_output_path], filtered_output_path)

    def test_Nick_testcase_2(self):
        """
        Test case using Nick's custom made maf and bam files for fillout testing

        This test cases uses the germline filter to exclude some mutations in the output

        This test case uses only research samples

        Takes about 10min to run
        """
        self.maxDiff = None
        self.runner_args['use_cache'] = False # do not use cache because it breaks for some reason
        self.runner_args['debug'] = True
        self.runner_args['js_console'] = True

        sample1_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample1.FillOutUnitTest01.muts.maf')
        sample2_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample2.FillOutUnitTest01.muts.maf')
        sample3_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample3.FillOutUnitTest01.muts.maf')
        sample4_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample4.FillOutUnitTest01.muts.maf')
        sample5_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample5.FillOutUnitTest01.muts.maf')

        sample1_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample1.UnitTest01.bam')
        sample2_bam =os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample2.UnitTest01.bam')
        sample3_bam =os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample3.UnitTest01.bam')
        sample4_bam =os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample4.UnitTest01.bam')
        sample5_bam =os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample5.UnitTest01.bam')

        self.input = {
            "samples": [
                {
                    "sample_id": "Sample1",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample1_maf },
                    "bam_file": { "class": "File", "path": sample1_bam }
                },
                {
                    "sample_id": "Sample2",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample2_maf },
                    "bam_file": { "class": "File", "path": sample2_bam }
                },
                {
                    "sample_id": "Sample3",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample3_maf },
                    "bam_file": { "class": "File", "path": sample3_bam }
                },
                {
                    "sample_id": "Sample4",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample4_maf },
                    "bam_file": { "class": "File", "path": sample4_bam }
                },
                {
                    "sample_id": "Sample5",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample5_maf },
                    "bam_file": { "class": "File", "path": sample5_bam }
                },
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']}
        }
        output_json, output_dir = self.run_cwl()
        output_path = os.path.join(output_dir,'output.maf')
        filtered_output_path = os.path.join(output_dir,'output.filtered.maf')
        portal_output_path = os.path.join(output_dir,'data_mutations_extended.txt')
        uncalled_output_path = os.path.join(output_dir,'data_mutations_uncalled.txt')

        expected_output = {
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        }

        # file contenst are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)
        # all_effects field is variable and changes bytes and checksum
        # need to check number of variant outputs instead
        self.assertNumMutationsHash(output_path, 475, 'd041bc641d85761b60c6b7ef8606bab2')
        self.assertNumMutationsHash(filtered_output_path, 475, 'd041bc641d85761b60c6b7ef8606bab2')
        self.assertNumMutationsHash(portal_output_path, 408, '63969ef90cb7a4524ab9063b4889bbde')
        self.assertNumMutationsHash(uncalled_output_path, 67, 'a474b61268d2a4c25fd27cc2ccbbce96')
        self.assertEqualNumMutations([portal_output_path, uncalled_output_path], filtered_output_path)

if __name__ == "__main__":
    unittest.main()
