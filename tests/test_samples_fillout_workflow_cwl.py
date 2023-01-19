#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the samples_fillout_workflow cwl
"""
import os
import sys
from datasets import (
    DATA_SETS,
)
from fixtures_fillout import rows
from pluto import (
    PlutoTestCase,
    PlutoPreRunTestCase,
    CWLFile,
    TableReader,
    md5_obj,
    ENABLE_LARGE_TESTS,
    OFile
)

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


class TestSamplesFilloutMixedGroup1(PlutoPreRunTestCase):
    cwl_file = CWLFile('samples_fillout_workflow.cwl')
    def setUpRun(self):
        """
        Run the workflow and return the results; output accessible under self.res.output in downstream 'test_' methods
        """
        sample_group1 = [
            {
                # THIS ONE LACKS A MAF FILE
                "sample_id": "Sample1",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                # "maf_file": { "class": "File", "path": sample1_maf },
                "bam_file": { "class": "File", "path": sample1_bam }
            },
            {
                "sample_id": "Sample2",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": sample2_maf }, # 66 muts
                "bam_file": { "class": "File", "path": sample2_bam }
            },
            {
                "sample_id": "Sample3",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": sample3_maf }, # 32 muts
                "bam_file": { "class": "File", "path": sample3_bam }
            },
        ]

        self.input = {
            "samples": sample_group1,
            "fillout_output_fname": 'output.maf',
            "ref_fasta": {"class": "File", "path": DATA_SETS['Proj_08390_G']['REF_FASTA']},
        }

        output_json, output_dir = self.run_cwl()

        return(output_json, output_dir)

    def getExpected(self, output_dir):
        """
        Return the expected CWL workflow output with the tmpdir output dir path included
        Accessible in downstream 'test_' methods under self.res.expected
        """
        return({
            "merged_vcf": OFile(name = "merged.vcf", dir = output_dir),
            "fillout_sources_vcf": OFile(name = "fillout.merged.sources.vcf", dir = output_dir),
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        })

    def test_CWLDictEqual(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        """
        # file contents are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', "fillout.merged.sources.vcf", ['size', 'checksum']),
        ('basename', 'merged.vcf', ['size', 'checksum']),
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(
                self.res.output,
                self.res.expected,
                related_keys = strip_related_keys)

    def test_output_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['output_file']).path, 258)

    def test_output_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['output_file']).path, "3a0f1991e37ab45885e8d5358f5a7ca6")


    def test_filtered_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['filtered_file']).path, 111)

    def test_filtered_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['filtered_file']).path, "a60d805c0d801d9f18174c5c5c96f841")

    def test_portal_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['portal_file']).path, 84)

    def test_portal_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['portal_file']).path, "742c992d5c04e1f3fda0d82eabf5f407")

    def test_uncalled_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['uncalled_file']).path, 27)

    def test_uncalled_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['uncalled_file']).path, "6ddbdc80180210ac4ad79050c72cfd89")

    def test_portal_output_path_num_muts(self):
        self.assertEqualNumMutations([
            OFile.init_dict(self.res.output['portal_file']).path,
            OFile.init_dict(self.res.output['uncalled_file']).path,
            ],
            OFile.init_dict(self.res.output['filtered_file']).path)

    def test_output_file_fields(self):
        self.assertMutFieldContains(
            OFile.init_dict(self.res.output['output_file']).path,
            "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3"], containsAll = True)

    def test_portal_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['portal_file']).path,
            "Amino_Acid_Change", [""])

    def test_uncalled_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['uncalled_file']).path,
            "Amino_Acid_Change", [""])





class TestSamplesFilloutMixedGroup2(PlutoPreRunTestCase):
    """
    Test that it runs when two samples are provided and one is lacking a maf because that triggers some of the singleton handling logic that we need to make sure works
    This test is super important!!
    """
    cwl_file = CWLFile('samples_fillout_workflow.cwl')
    def setUpRun(self):
        """
        Run the workflow and return the results; output accessible under self.res.output in downstream 'test_' methods
        """
        sample_group1 = [
            {
                # THIS ONE LACKS A MAF FILE
                "sample_id": "Sample1",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                # "maf_file": { "class": "File", "path": sample1_maf },
                "bam_file": { "class": "File", "path": sample1_bam }
            },
            {
                "sample_id": "Sample2",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": sample2_maf },
                "bam_file": { "class": "File", "path": sample2_bam }
            },
        ]

        self.input = {
            "samples": sample_group1,
            "fillout_output_fname": 'output.maf',
            "ref_fasta": {"class": "File", "path": DATA_SETS['Proj_08390_G']['REF_FASTA']},
        }

        output_json, output_dir = self.run_cwl()

        return(output_json, output_dir)

    def getExpected(self, output_dir):
        """
        Return the expected CWL workflow output with the tmpdir output dir path included
        Accessible in downstream 'test_' methods under self.res.expected
        """
        return({
            "merged_vcf": OFile(name = "Sample2.sorted.vcf", dir = output_dir),
            "fillout_sources_vcf": OFile(name = "fillout.merged.sources.vcf", dir = output_dir),
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        })

    def test_CWLDictEqual(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        """
        # file contents are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', "fillout.merged.sources.vcf", ['size', 'checksum']),
        ('basename', 'Sample2.sorted.vcf', ['size', 'checksum']),
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(
                self.res.output,
                self.res.expected,
                related_keys = strip_related_keys)

    def test_output_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['output_file']).path, 132)

    def test_output_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['output_file']).path, "1cbf98a55bb79e2b304ed8966e7282d1")


    def test_filtered_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['filtered_file']).path, 132)

    def test_filtered_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['filtered_file']).path, "1cbf98a55bb79e2b304ed8966e7282d1")

    def test_portal_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['portal_file']).path, 132)

    def test_portal_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['portal_file']).path, "c54e74773a7203e4823ca47ed344e4e6")

    def test_uncalled_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['uncalled_file']).path, 0)

    def test_uncalled_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['uncalled_file']).path, "d751713988987e9331980363e24189ce")

    def test_portal_output_path_num_muts(self):
        self.assertEqualNumMutations([
            OFile.init_dict(self.res.output['portal_file']).path,
            OFile.init_dict(self.res.output['uncalled_file']).path,
            ],
            OFile.init_dict(self.res.output['filtered_file']).path)

    def test_output_file_fields(self):
        self.assertMutFieldContains(
            OFile.init_dict(self.res.output['output_file']).path,
            "Tumor_Sample_Barcode", ["Sample1", "Sample2"], containsAll = True)

    def test_portal_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['portal_file']).path,
            "Amino_Acid_Change", [""])

    def test_uncalled_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['uncalled_file']).path,
            "Amino_Acid_Change", [""])







class TestSamplesFillout1(PlutoPreRunTestCase):
    """
    test_Nick_testcase
    Test case using Nick's custom made maf and bam files for fillout testing

    This test cases uses the germline filter to exclude some mutations in the output

    Takes about 10min to run
    """
    cwl_file = CWLFile('samples_fillout_workflow.cwl')
    def setUpRun(self):
        """
        Run the workflow and return the results; output accessible under self.res.output in downstream 'test_' methods
        """
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
            "ref_fasta": {"class": "File", "path": DATA_SETS['Proj_08390_G']['REF_FASTA']}
        }

        output_json, output_dir = self.run_cwl()

        return(output_json, output_dir)

    def getExpected(self, output_dir):
        """
        Return the expected CWL workflow output with the tmpdir output dir path included
        Accessible in downstream 'test_' methods under self.res.expected
        """
        return({
            "merged_vcf": OFile(name = "merged.vcf", dir = output_dir),
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            "fillout_sources_vcf": OFile(name = "fillout.merged.sources.vcf", dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        })

    def test_CWLDictEqual(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        """
        # file contents are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', "fillout.merged.sources.vcf", ['size', 'checksum']),
        ('basename', 'merged.vcf', ['size', 'checksum']),
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(
                self.res.output,
                self.res.expected,
                related_keys = strip_related_keys)

    def test_output_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['output_file']).path, 475)

    def test_output_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['output_file']).path, "d041bc641d85761b60c6b7ef8606bab2")


    def test_filtered_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['filtered_file']).path, 230)

    def test_filtered_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['filtered_file']).path, "c9cde01507d1b2470057c5d120eaab68")

    def test_portal_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['portal_file']).path, 163)

    def test_portal_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['portal_file']).path, "2bdb5afc25d0f4af2fb53cb9edfdc03f")

    def test_uncalled_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['uncalled_file']).path, 67)

    def test_uncalled_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['uncalled_file']).path, "42f32ce8b3c4ffff4a8d26bd9e3bb900")

    def test_portal_output_path_num_muts(self):
        self.assertEqualNumMutations([
            OFile.init_dict(self.res.output['portal_file']).path,
            OFile.init_dict(self.res.output['uncalled_file']).path,
            ],
            OFile.init_dict(self.res.output['filtered_file']).path)

    def test_output_file_fields(self):
        self.assertMutFieldContains(
            OFile.init_dict(self.res.output['output_file']).path,
            "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3", 'Sample4', 'Sample5'], containsAll = True)

    def test_portal_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['portal_file']).path,
            "Amino_Acid_Change", [""])

    def test_uncalled_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['uncalled_file']).path,
            "Amino_Acid_Change", [""])




class TestSamplesFillout2(PlutoPreRunTestCase):
    """
    test_Nick_testcase_2
    Test case using Nick's custom made maf and bam files for fillout testing

    This test cases uses the germline filter to exclude some mutations in the output

    This test case uses only research samples

    Takes about 10min to run

    """
    cwl_file = CWLFile('samples_fillout_workflow.cwl')
    cwl_file = CWLFile('samples_fillout_workflow.cwl')
    def setUpRun(self):
        """
        Run the workflow and return the results; output accessible under self.res.output in downstream 'test_' methods
        """
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
            "ref_fasta": {"class": "File", "path": DATA_SETS['Proj_08390_G']['REF_FASTA']}
        }

        output_json, output_dir = self.run_cwl()

        return(output_json, output_dir)

    def getExpected(self, output_dir):
        """
        Return the expected CWL workflow output with the tmpdir output dir path included
        Accessible in downstream 'test_' methods under self.res.expected
        """
        return({
            "merged_vcf": OFile(name = "merged.vcf", dir = output_dir),
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            "fillout_sources_vcf": OFile(name = "fillout.merged.sources.vcf", dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        })

    def test_CWLDictEqual(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        """
        # file contents are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', "fillout.merged.sources.vcf", ['size', 'checksum']),
        ('basename', 'merged.vcf', ['size', 'checksum']),
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(
                self.res.output,
                self.res.expected,
                related_keys = strip_related_keys)

    def test_output_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['output_file']).path, 475)

    def test_output_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['output_file']).path, "d041bc641d85761b60c6b7ef8606bab2")


    def test_filtered_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['filtered_file']).path, 475)

    def test_filtered_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['filtered_file']).path, "d041bc641d85761b60c6b7ef8606bab2")

    def test_portal_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['portal_file']).path, 408)

    def test_portal_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['portal_file']).path, "0bcea8598de4d5816aa44aa63c6790ca")

    def test_uncalled_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['uncalled_file']).path, 67)

    def test_uncalled_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['uncalled_file']).path, "42f32ce8b3c4ffff4a8d26bd9e3bb900")

    def test_portal_output_path_num_muts(self):
        self.assertEqualNumMutations([
            OFile.init_dict(self.res.output['portal_file']).path,
            OFile.init_dict(self.res.output['uncalled_file']).path,
            ],
            OFile.init_dict(self.res.output['filtered_file']).path)

    def test_output_file_fields(self):
        self.assertMutFieldContains(
            OFile.init_dict(self.res.output['output_file']).path,
            "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3", 'Sample4', 'Sample5'], containsAll = True)

    def test_portal_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['portal_file']).path,
            "Amino_Acid_Change", [""])

    def test_uncalled_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['uncalled_file']).path,
            "Amino_Acid_Change", [""])
