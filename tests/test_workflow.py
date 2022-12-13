#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the workflow.cwl
"""
import os
import sys
from datasets import (
    DATA_SETS,
    KNOWN_FUSIONS_FILE,
    IMPACT_FILE,
)
from pluto import (
    PlutoTestCase,
    OFile,
    ODir,
    serialize_repr,
    load_mutations,
    run_cwl,
    CWLFile,
)



class TestWorkflow(PlutoTestCase):
    cwl_file = CWLFile('workflow.cwl')

    def test_run_worflow_one_maf(self):
        """
        Test that the workflow works correctly when run with a single maf
        """
        self.skipTest("Fix jenkins error")
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        self.input = {
            "project_id": "Proj_08390_G",
            "project_name": "Proj_08390_G",
            "project_short_name": "Proj_08390_G",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_08390_G',
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "analysis_mutations_share_filename": "Proj_08390_G.muts.share.maf",
            "analysis_segment_cna_filename": "Proj_08390_G.seg.cna.txt",
            "analysis_sv_filename": "Proj_08390_G.svs.maf",
            "cbio_meta_cna_segments_filename": "Proj_08390_G_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_08390_G_data_cna_hg19.seg",
            "helix_filter_version": "20.06.1",
            'IMPACT_gene_list': {
                "path": IMPACT_FILE,
                "class": "File"
            },
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
            },
            "sample_summary_file": {
                "path": sample_summary_file,
                "class": "File"
            },
            "targets_list": {
                "path": DATA_SETS['Proj_08390_G']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                    "class": "File"
                }
            ],
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "mutation_svs_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.maf"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                }
            ],
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
        'analysis_dir': ODir(name='analysis', items=[
            OFile(name='Proj_08390_G.gene.cna.txt', size=87905, hash='7cc89d24556de93b9a409812317581e67e5df494'),
            OFile(name='Proj_08390_G.muts.maf', size=33243, hash='2c8904927a917d6e935ef207582d995680574d16'),
            OFile(name='Proj_08390_G.muts.share.maf', size=7462, hash='b5af4e0fcd89fecabf8095aa3d7690e5edb8dca1'),
            OFile(name='Proj_08390_G.seg.cna.txt', size=1632, hash='f0ebb82c34b6530447fa1e70b6dedcc039840d61'),
            OFile(name='Proj_08390_G.svs.maf', size=23603, hash='df420706bb5b772a79317843c0a01a3c88a9571d')], dir=output_dir),
        'portal_dir': ODir(name='portal', items=[
            OFile(name='meta_clinical_sample.txt', size=140, hash='4c567d81c3b17a76c324fd3e2f73793a6e804f65'),
            OFile(name='data_clinical_patient.txt', size=643, hash='9417dcabddd6ab2cbe98167bccd9b9e4fa182562'),
            OFile(name='data_clinical_sample.txt', size=7592, hash='2a0c59593fa7726743b2fe46db9d955dbc625453'),
            OFile(name='meta_study.txt', size=152, hash='2b0a5fd1a97329adf7c3b1596c84cd6567059a95'),
            OFile(name='meta_clinical_patient.txt', size=142, hash='9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f'),
            OFile(name='meta_CNA.txt', size=270, hash='a9bf16f6a0490b19e611e8814b85f7bf1d52417a'),
            OFile(name='meta_fusions.txt', size=227, hash='77649e888bafc6a4ed61261d1c46d2f238e1c32b'),
            OFile(name='meta_mutations_extended.txt', size=253, hash='fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d'),
            OFile(name='Proj_08390_G_meta_cna_hg19_seg.txt', size=200, hash='59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e'),
            OFile(name='data_CNA.txt', size=5365, hash='931d82412733d7f93dd4117cd955f35e5dcbacc1'),
            OFile(name='data_CNA.ascna.txt', size=6164, hash='452d5ddef12a44693d5a98a05f5d300801734cfe'),
            OFile(name='data_mutations_extended.txt', size=5106, hash='e713516cf04750a3e3f1ef932b1c7202d4b75bf2'),
            OFile(name='Proj_08390_G_data_cna_hg19.seg', size=1632, hash='f0ebb82c34b6530447fa1e70b6dedcc039840d61'),
            OFile(name='data_fusions.txt', size=99, hash='c16f763b248813fcdde76f7486f1ddc4e9856038'),
            ODir(name='case_lists', items=[
                OFile(name='cases_all.txt', size=616, hash='b9e43289cec5603b0886b5e8507c8d019387c125'),
                OFile(name='cases_cnaseq.txt', size=696, hash='b87e2da8dce0fddbadec348efe2986519b2a794b'),
                OFile(name='cases_cna.txt', size=628, hash='053481a8299e9430117f8e45e081aa7ec21033a6'),
                OFile(name='cases_sequenced.txt', size=641, hash='ef9f5aef03c2527bf576470168660557ca1c7cc9')])],
                dir=output_dir)}

        self.assertCWLDictEqual(output_json, expected_output)
        comments, mutations = load_mutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'))
        self.assertEqual(len(mutations), 22)
        comments, mutations = load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 17)

        # load the data_CNA.txt file
        path = os.path.join(output_dir, 'portal/data_CNA.txt') # renamed from the data_CNA.scna.txt file ...
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 's_C_VJ7F47_P001_d']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 's_C_VJ7F47_P001_d']
        self.assertEqual(header_parts, expected_header_parts)



    def test_run_worflow_two_mafs(self):
        """
        Test that the workflow works correctly when run with two maf files
        """
        self.skipTest("Fix jenkins error")
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        self.input = {
            "project_id": "Proj_08390_G",
            "project_name": "Proj_08390_G",
            "project_short_name": "Proj_08390_G",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_08390_G',
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "analysis_mutations_share_filename": "Proj_08390_G.muts.share.maf",
            "analysis_segment_cna_filename": "Proj_08390_G.seg.cna.txt",
            "analysis_sv_filename": "Proj_08390_G.svs.maf",
            "cbio_meta_cna_segments_filename": "Proj_08390_G_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_08390_G_data_cna_hg19.seg",
            "helix_filter_version": "20.06.1",
            'IMPACT_gene_list': {
                "path": IMPACT_FILE,
                "class": "File"
            },
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
            },
            "sample_summary_file": {
                "path": sample_summary_file,
                "class": "File"
            },
            "targets_list": {
                "path": DATA_SETS['Proj_08390_G']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.muts.maf"),
                    "class": "File"
                }
            ],
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "mutation_svs_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.maf"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.svs.pass.vep.maf"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample3.rg.md.abra.printreads__Sample4.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample3.rg.md.abra.printreads__Sample4.rg.md.abra.printreads_hisens.seg"),
                    "class": "File"
                }
            ]
        }

        output_json, output_dir = self.run_cwl()

        self.maxDiff = None
        expected_output = {
        'analysis_dir': ODir(name='analysis', items=[
            OFile(name='Proj_08390_G.gene.cna.txt', size=173982, hash='ab17d587ad5ae0a87fd6c6d4dd2d5d1701208ce9'),
            OFile(name='Proj_08390_G.muts.maf', size=54458, hash='d4352ee2b702877b84db2b632972ccad2441f3e0'),
            OFile(name='Proj_08390_G.muts.share.maf', size=10956, hash='086ce6517eae68e47160c8740c5f00d7c3454110'),
            OFile(name='Proj_08390_G.seg.cna.txt', size=3191, hash='f6a77b280c047a7e2082e3a09e8138f861790d3a'),
            OFile(name='Proj_08390_G.svs.maf', size=35595, hash='5c2a63fc01980550108e58079a8b689d53c97d8c')], dir=output_dir),
        'portal_dir': ODir(name='portal', items=[
            OFile(name='meta_clinical_sample.txt', size=140, hash='4c567d81c3b17a76c324fd3e2f73793a6e804f65'),
            OFile(name='data_clinical_patient.txt', size=643, hash='9417dcabddd6ab2cbe98167bccd9b9e4fa182562'),
            OFile(name='data_clinical_sample.txt', size=7592, hash='2a0c59593fa7726743b2fe46db9d955dbc625453'),
            OFile(name='meta_study.txt', size=152, hash='2b0a5fd1a97329adf7c3b1596c84cd6567059a95'),
            OFile(name='meta_clinical_patient.txt', size=142, hash='9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f'),
            OFile(name='meta_CNA.txt', size=270, hash='a9bf16f6a0490b19e611e8814b85f7bf1d52417a'),
            OFile(name='meta_fusions.txt', size=227, hash='77649e888bafc6a4ed61261d1c46d2f238e1c32b'),
            OFile(name='meta_mutations_extended.txt', size=253, hash='fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d'),
            OFile(name='Proj_08390_G_meta_cna_hg19_seg.txt', size=200, hash='59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e'),
            OFile(name='data_CNA.txt', size=6784, hash='09b4d944e50ea9d0e7567e04ce55b0f21d281255'),
            OFile(name='data_CNA.ascna.txt', size=8789, hash='d93ffe83137d9a77e2420b40ab3a2e0a1a5ad069'),
            OFile(name='data_mutations_extended.txt', size=7539, hash='43469aa0f9125d3dca6217ee02641638c3a92e24'),
            OFile(name='Proj_08390_G_data_cna_hg19.seg', size=3191, hash='f6a77b280c047a7e2082e3a09e8138f861790d3a'),
            OFile(name='data_fusions.txt', size=99, hash='c16f763b248813fcdde76f7486f1ddc4e9856038'),
            ODir(name='case_lists', items=[
                OFile(name='cases_all.txt', size=616, hash='b9e43289cec5603b0886b5e8507c8d019387c125'),
                OFile(name='cases_cnaseq.txt', size=696, hash='b87e2da8dce0fddbadec348efe2986519b2a794b'),
                OFile(name='cases_cna.txt', size=628, hash='053481a8299e9430117f8e45e081aa7ec21033a6'),
                OFile(name='cases_sequenced.txt', size=641, hash='ef9f5aef03c2527bf576470168660557ca1c7cc9')])],
                dir=output_dir)}
        self.assertCWLDictEqual(output_json, expected_output)

        comments, mutations = load_mutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'))
        self.assertEqual(len(mutations), 34)
        comments, mutations = load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 27)

        # TODO: fix the Facets column headers so that this passes
        # load the data_CNA.txt file
        path = os.path.join(output_dir, 'portal/data_CNA.txt') # renamed from the data_CNA.scna.txt file ...
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 's_C_VJ7F47_P001_d', 's_C_X50T9Y_P001_d']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 's_C_VJ7F47_P001_d', 's_C_X50T9Y_P001_d']
        self.assertEqual(header_parts, expected_header_parts)
