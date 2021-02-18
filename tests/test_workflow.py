#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the workflow.cwl
"""
import os
import sys
import unittest
from tempfile import TemporaryDirectory

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import load_mutations, run_cwl, CWLFile
from pluto.settings import DATA_SETS, KNOWN_FUSIONS_FILE, IMPACT_FILE
sys.path.pop(0)

cwl_file = CWLFile('workflow.cwl')

class TestWorkflow(unittest.TestCase):
    def test_run_worflow_one_maf(self):
        """
        Test that the workflow works correctly when run with a single maf
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        input_json = {
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

        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

            expected_output = {
                'analysis_dir': {
                    'class': 'Directory',
                    'basename': 'analysis',
                    'location': 'file://' + os.path.join(output_dir, 'analysis'),
                    'path': os.path.join(output_dir, 'analysis'),
                    'listing': [
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt'),
                            'basename': 'Proj_08390_G.gene.cna.txt',
                            'class': 'File',
                            'checksum': 'sha1$7cc89d24556de93b9a409812317581e67e5df494',
                            'size': 87905,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                            'basename': 'Proj_08390_G.muts.maf',
                            'class': 'File',
                            'checksum': 'sha1$2c8904927a917d6e935ef207582d995680574d16',
                            'size': 33243,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf'),
                            'basename': 'Proj_08390_G.muts.share.maf',
                            'class': 'File',
                            'checksum': 'sha1$b5af4e0fcd89fecabf8095aa3d7690e5edb8dca1',
                            'size': 7462,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt'),
                            'basename': 'Proj_08390_G.seg.cna.txt',
                            'class': 'File',
                            'checksum': 'sha1$f0ebb82c34b6530447fa1e70b6dedcc039840d61',
                            'size': 1632,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf'),
                            'basename': 'Proj_08390_G.svs.maf',
                            'class': 'File',
                            'checksum': 'sha1$df420706bb5b772a79317843c0a01a3c88a9571d',
                            'size': 23603,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf')
                            }
                        ]
                    },
                    'portal_dir': {
                        'location': 'file://' + os.path.join(output_dir, 'portal'),
                        'path': os.path.join(output_dir, 'portal'),
                        'class': 'Directory',
                        'basename': 'portal',
                        'listing': [
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_sample.txt'),
                                'basename': 'meta_clinical_sample.txt',
                                'class': 'File',
                                'checksum': 'sha1$4c567d81c3b17a76c324fd3e2f73793a6e804f65',
                                'size': 140,
                                'path': os.path.join(output_dir, 'portal/meta_clinical_sample.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_patient.txt'),
                                'basename': 'data_clinical_patient.txt',
                                'class': 'File',
                                'checksum': 'sha1$9417dcabddd6ab2cbe98167bccd9b9e4fa182562',
                                'size': 643,
                                'path': os.path.join(output_dir, 'portal/data_clinical_patient.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_sample.txt'),
                                'basename': 'data_clinical_sample.txt',
                                'class': 'File',
                                'checksum': 'sha1$2a0c59593fa7726743b2fe46db9d955dbc625453',
                                'size': 7592,
                                'path': os.path.join(output_dir, 'portal/data_clinical_sample.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_study.txt'),
                                'basename': 'meta_study.txt',
                                'class': 'File',
                                'checksum': 'sha1$2b0a5fd1a97329adf7c3b1596c84cd6567059a95',
                                'size': 152,
                                'path': os.path.join(output_dir, 'portal/meta_study.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_patient.txt'),
                                'basename': 'meta_clinical_patient.txt',
                                'class': 'File',
                                'checksum': 'sha1$9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f',
                                'size': 142,
                                'path': os.path.join(output_dir, 'portal/meta_clinical_patient.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_CNA.txt'),
                                'basename': 'meta_CNA.txt',
                                'class': 'File',
                                'checksum': 'sha1$a9bf16f6a0490b19e611e8814b85f7bf1d52417a',
                                'size': 270,
                                'path': os.path.join(output_dir, 'portal/meta_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_fusions.txt'),
                                'basename': 'meta_fusions.txt',
                                'class': 'File',
                                'checksum': 'sha1$77649e888bafc6a4ed61261d1c46d2f238e1c32b',
                                'size': 227,
                                'path': os.path.join(output_dir, 'portal/meta_fusions.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_mutations_extended.txt'),
                                'basename': 'meta_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d',
                                'size': 253,
                                'path': os.path.join(output_dir, 'portal/meta_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt'),
                                'basename': 'Proj_08390_G_meta_cna_hg19_seg.txt',
                                'class': 'File',
                                'checksum': 'sha1$59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e',
                                'size': 200,
                                'path': os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.txt'),
                                'basename': 'data_CNA.txt',
                                'class': 'File',
                                'checksum': 'sha1$931d82412733d7f93dd4117cd955f35e5dcbacc1',
                                'size': 5365,
                                'path': os.path.join(output_dir, 'portal/data_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                                'basename': 'data_CNA.ascna.txt',
                                'class': 'File',
                                'checksum': 'sha1$452d5ddef12a44693d5a98a05f5d300801734cfe',
                                'size': 6164,
                                'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                                'basename': 'data_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$e713516cf04750a3e3f1ef932b1c7202d4b75bf2',
                                'size': 5106,
                                'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                                'basename': 'Proj_08390_G_data_cna_hg19.seg',
                                'class': 'File',
                                'checksum': 'sha1$f0ebb82c34b6530447fa1e70b6dedcc039840d61',
                                'size': 1632,
                                'path': os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_fusions.txt'),
                                'basename': 'data_fusions.txt',
                                'class': 'File',
                                'checksum': 'sha1$c16f763b248813fcdde76f7486f1ddc4e9856038',
                                'size': 99,
                                'path': os.path.join(output_dir, 'portal/data_fusions.txt')
                            },
                            {
                                'class': 'Directory',
                                'basename': 'case_lists',
                                'location': 'file://' + os.path.join(output_dir, 'portal/case_lists'),
                                'path': os.path.join(output_dir, 'portal/case_lists'),
                                'listing': [
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_all.txt'),
                                    'basename': 'cases_all.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$b9e43289cec5603b0886b5e8507c8d019387c125',
                                    'size': 616,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_all.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt'),
                                    'basename': 'cases_cnaseq.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$b87e2da8dce0fddbadec348efe2986519b2a794b',
                                    'size': 696,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cna.txt'),
                                    'basename': 'cases_cna.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$053481a8299e9430117f8e45e081aa7ec21033a6',
                                    'size': 628,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_cna.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt'),
                                    'basename': 'cases_sequenced.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$ef9f5aef03c2527bf576470168660557ca1c7cc9',
                                    'size': 641,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt')}
                                ]
                            }
                        ]
                    }
                }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)
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
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        input_json = {
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
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

            expected_output = {
                'analysis_dir': {
                    'class': 'Directory',
                    'basename': 'analysis',
                    'location': 'file://' + os.path.join(output_dir, 'analysis'),
                    'path': os.path.join(output_dir, 'analysis'),
                    'listing': [
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt'),
                        'basename': 'Proj_08390_G.gene.cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$ab17d587ad5ae0a87fd6c6d4dd2d5d1701208ce9',
                        'size': 173982,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                        'basename': 'Proj_08390_G.muts.maf',
                        'class': 'File',
                        'checksum': 'sha1$d4352ee2b702877b84db2b632972ccad2441f3e0',
                        'size': 54458,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')},
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf'),
                            'basename': 'Proj_08390_G.muts.share.maf',
                            'class': 'File',
                            'checksum': 'sha1$086ce6517eae68e47160c8740c5f00d7c3454110',
                            'size': 10956,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf')
                        },
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt'),
                        'basename': 'Proj_08390_G.seg.cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$f6a77b280c047a7e2082e3a09e8138f861790d3a',
                        'size': 3191,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf'),
                        'basename': 'Proj_08390_G.svs.maf',
                        'class': 'File',
                        'checksum': 'sha1$5c2a63fc01980550108e58079a8b689d53c97d8c',
                        'size': 35595,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf')}
                    ]
                },
                'portal_dir': {
                    'class': 'Directory',
                    'basename': 'portal',
                    'location': 'file://' + os.path.join(output_dir, 'portal'),
                    'path': os.path.join(output_dir, 'portal'),
                    'listing': [
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_sample.txt'),
                        'basename': 'meta_clinical_sample.txt',
                        'class': 'File',
                        'checksum': 'sha1$4c567d81c3b17a76c324fd3e2f73793a6e804f65',
                        'size': 140,
                        'path': os.path.join(output_dir, 'portal/meta_clinical_sample.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_patient.txt'),
                        'basename': 'data_clinical_patient.txt',
                        'class': 'File',
                        'checksum': 'sha1$9417dcabddd6ab2cbe98167bccd9b9e4fa182562',
                        'size': 643,
                        'path': os.path.join(output_dir, 'portal/data_clinical_patient.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_sample.txt'),
                        'basename': 'data_clinical_sample.txt',
                        'class': 'File',
                        'checksum': 'sha1$2a0c59593fa7726743b2fe46db9d955dbc625453',
                        'size': 7592,
                        'path': os.path.join(output_dir, 'portal/data_clinical_sample.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_study.txt'),
                        'basename': 'meta_study.txt',
                        'class': 'File',
                        'checksum': 'sha1$2b0a5fd1a97329adf7c3b1596c84cd6567059a95',
                        'size': 152,
                        'path': os.path.join(output_dir, 'portal/meta_study.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_patient.txt'),
                        'basename': 'meta_clinical_patient.txt',
                        'class': 'File',
                        'checksum': 'sha1$9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f',
                        'size': 142,
                        'path': os.path.join(output_dir, 'portal/meta_clinical_patient.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_CNA.txt'),
                        'basename': 'meta_CNA.txt',
                        'class': 'File',
                        'checksum': 'sha1$a9bf16f6a0490b19e611e8814b85f7bf1d52417a',
                        'size': 270,
                        'path': os.path.join(output_dir, 'portal/meta_CNA.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_fusions.txt'),
                        'basename': 'meta_fusions.txt',
                        'class': 'File',
                        'checksum': 'sha1$77649e888bafc6a4ed61261d1c46d2f238e1c32b',
                        'size': 227,
                        'path': os.path.join(output_dir, 'portal/meta_fusions.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_mutations_extended.txt'),
                        'basename': 'meta_mutations_extended.txt',
                        'class': 'File',
                        'checksum': 'sha1$fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d',
                        'size': 253,
                        'path': os.path.join(output_dir, 'portal/meta_mutations_extended.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt'),
                        'basename': 'Proj_08390_G_meta_cna_hg19_seg.txt',
                        'class': 'File',
                        'checksum': 'sha1$59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e',
                        'size': 200,
                        'path': os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.txt'),
                        'basename': 'data_CNA.txt',
                        'class': 'File',
                        'checksum': 'sha1$09b4d944e50ea9d0e7567e04ce55b0f21d281255',
                        'size': 6784,
                        'path': os.path.join(output_dir, 'portal/data_CNA.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                        'basename': 'data_CNA.ascna.txt',
                        'class': 'File',
                        'checksum': 'sha1$d93ffe83137d9a77e2420b40ab3a2e0a1a5ad069',
                        'size': 8789,
                        'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                        'basename': 'data_mutations_extended.txt',
                        'class': 'File',
                        'checksum': 'sha1$43469aa0f9125d3dca6217ee02641638c3a92e24',
                        'size': 7539,
                        'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                        'basename': 'Proj_08390_G_data_cna_hg19.seg',
                        'class': 'File',
                        'checksum': 'sha1$f6a77b280c047a7e2082e3a09e8138f861790d3a',
                        'size': 3191,
                        'path': os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_fusions.txt'),
                        'basename': 'data_fusions.txt',
                        'class': 'File',
                        'checksum': 'sha1$c16f763b248813fcdde76f7486f1ddc4e9856038',
                        'size': 99,
                        'path': os.path.join(output_dir, 'portal/data_fusions.txt')},
                        {'class': 'Directory',
                        'basename': 'case_lists',
                        'listing': [{'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_all.txt'),
                        'basename': 'cases_all.txt',
                        'class': 'File',
                        'checksum': 'sha1$b9e43289cec5603b0886b5e8507c8d019387c125',
                        'size': 616,
                        'path': os.path.join(output_dir, 'portal/case_lists/cases_all.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt'),
                        'basename': 'cases_cnaseq.txt',
                        'class': 'File',
                        'checksum': 'sha1$b87e2da8dce0fddbadec348efe2986519b2a794b',
                        'size': 696,
                        'path': os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cna.txt'),
                        'basename': 'cases_cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$053481a8299e9430117f8e45e081aa7ec21033a6',
                        'size': 628,
                        'path': os.path.join(output_dir, 'portal/case_lists/cases_cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt'),
                        'basename': 'cases_sequenced.txt',
                        'class': 'File',
                        'checksum': 'sha1$ef9f5aef03c2527bf576470168660557ca1c7cc9',
                        'size': 641,
                        'path': os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt')}],
                        'location': 'file://' + os.path.join(output_dir, 'portal/case_lists'),
                        'path': os.path.join(output_dir, 'portal/case_lists')}
                    ]
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)
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

if __name__ == "__main__":
    unittest.main()
