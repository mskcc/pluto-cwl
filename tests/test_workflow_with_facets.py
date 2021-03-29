#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the workflow_with_facets.cwl
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import TableReader, PlutoTestCase
sys.path.pop(0)

class TestWorkflowWithFacets(PlutoTestCase):
    cwl_file = 'workflow_with_facets.cwl'

    def test_demo_dataset1(self):
        """
        Test case for using a single demo sample; tiny dataset
        """
        # self.preserve = True
        # print(self.tmpdir)
        data_clinical_file = os.path.join(self.DATA_SETS['demo']['INPUTS_DIR'], "demo_sample_data_clinical.txt")
        sample_summary_file = os.path.join(self.DATA_SETS['demo']['QC_DIR'], "demo_SampleSummary.txt")
        mutation_svs_txt_file = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt")
        mutation_svs_maf = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.maf")
        pair_maf = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        snp_pileup = os.path.join(self.DATA_SETS['demo']['SNP_PILEUP_DIR'], "Sample1.Sample2.snp_pileup.gz")

        self.input = {
            "assay_coverage": "10000", # TODO: get this from an assay reference key
            "project_id": "demo",
            "project_name": "demo",
            "project_short_name": "demo",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'demo',
            "analysis_gene_cna_filename": "demo.gene.cna.txt",
            "analysis_mutations_filename": "demo.muts.maf",
            "analysis_mutations_share_filename": "demo.muts.share.maf",
            "analysis_segment_cna_filename": "demo.seg.cna.txt",
            "analysis_sv_filename": "demo.svs.maf",
            "cbio_meta_cna_segments_filename": "Proj_08390_G_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_08390_G_data_cna_hg19.seg",
            "helix_filter_version": "20.06.1",
            'IMPACT_gene_list': {
                "path": self.IMPACT_FILE,
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
                "path": self.DATA_SETS['demo']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": self.KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_svs_txt_files": [ { "path": mutation_svs_txt_file, "class": "File" } ],
            "mutation_svs_maf_files": [ { "path": mutation_svs_maf, "class": "File" } ],
            "pairs": [
                {
                    "pair_maf": { "path": pair_maf, "class": "File" },
                    "snp_pileup": { "path": snp_pileup, "class": "File" },
                    "pair_id": "Sample1.Sample2",
                    "tumor_id": "Sample1",
                    "normal_id": "Sample2"
                }
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_dir': {
                'class': 'Directory',
                'basename': 'analysis',
                'location': 'file://' + os.path.join(output_dir, 'analysis'),
                'path': os.path.join(output_dir, 'analysis'),
                'listing': [
                    {
                        'location': 'file://' + os.path.join(output_dir, 'analysis/demo.gene.cna.txt'),
                        'basename': 'demo.gene.cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$a3527d2bf248c3121f093f7edfdf2605b12dfe71',
                        'size': 26505,
                        'path': os.path.join(output_dir, 'analysis/demo.gene.cna.txt')
                    },
                    {
                        'location': 'file://' + os.path.join(output_dir, 'analysis/demo.muts.maf'),
                        'basename': 'demo.muts.maf',
                        'class': 'File',
                        'checksum': 'sha1$bbb1320d792c5a37f0d2afae82fdd3e67b2db182',
                        'size': 36333,
                        'path': os.path.join(output_dir, 'analysis/demo.muts.maf')
                    },
                    {
                        'location': 'file://' + os.path.join(output_dir, 'analysis/demo.muts.share.maf'),
                        'basename': 'demo.muts.share.maf',
                        'class': 'File',
                        'checksum': 'sha1$7716a2274203abe8039c99d0722b92a675ed38b7',
                        'size': 7288,
                        'path': os.path.join(output_dir, 'analysis/demo.muts.share.maf')
                    },
                    {
                        'location': 'file://' + os.path.join(output_dir, 'analysis/demo.seg.cna.txt'),
                        'basename': 'demo.seg.cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$e99e98d33598c423448cba44a37da712c9c933c5',
                        'size': 525,
                        'path': os.path.join(output_dir, 'analysis/demo.seg.cna.txt')
                    },
                    {
                        'location': 'file://' + os.path.join(output_dir, 'analysis/demo.svs.maf'),
                        'basename': 'demo.svs.maf',
                        'class': 'File',
                        'checksum': 'sha1$df420706bb5b772a79317843c0a01a3c88a9571d',
                        'size': 23603,
                        'path': os.path.join(output_dir, 'analysis/demo.svs.maf')
                        }
                    ]
                },
                'facets_dir': {
                    'location': 'file://' + os.path.join(output_dir, 'facets'),
                    'basename': 'facets',
                    'class': 'Directory',
                    'path': os.path.join(output_dir, 'facets'),
                    'listing': [{'basename': 'Sample1.Sample2',
                          'class': 'Directory',
                          'listing': [{'basename': 'Sample1.Sample2_hisens.ccf.portal.maf',
                                       'checksum': 'sha1$e86566d281bc808e57ee912d02bdbda8179d7c0a',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.Sample2_hisens.ccf.portal.maf'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.Sample2_hisens.ccf.portal.maf'),
                                       'size': 62216},
                                      {'basename': 'Sample1.arm_level.txt',
                                       'checksum': 'sha1$98704689aedd071919ca0b2e387462948cb803b4',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.arm_level.txt'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.arm_level.txt'),
                                       'size': 370},
                                      {'basename': 'Sample1.txt',
                                       'checksum': 'sha1$ccdaded463980977bc4c90e539043c6efd1aef33',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.txt'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.txt'),
                                       'size': 506},
                                      {'basename': 'Sample1.gene_level.txt',
                                       'checksum': 'sha1$75ad4054ccec1f7a47d1c7738730816191d713c3',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.gene_level.txt'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.gene_level.txt'),
                                       'size': 94824},
                                      {'basename': 'Sample1_hisens.cncf.txt',
                                       'checksum': 'sha1$738974eaa3e4278f40cdbbf487b87f8b57ee2ecb',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.cncf.txt'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.cncf.txt'),
                                       'size': 1385},
                                      {'basename': 'Sample1_hisens.rds',
                                       'checksum': 'sha1$13846db5dc75ca6b7d64cb69ff04fa64da9a88e7',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.rds'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.rds'),
                                       'size': 191089},
                                      {'basename': 'Sample1_hisens.seg',
                                       'checksum': 'sha1$069b23dada099e20700c6b10471750e3cba6ee86',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.seg'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.seg'),
                                       'size': 647},
                                      {'basename': 'Sample1_purity.rds',
                                       'checksum': 'sha1$ded479ac6278f1cd7bf42bf186440d8a89415a02',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.rds'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.rds'),
                                       'size': 191195},
                                      {'basename': 'Sample1_purity.seg',
                                       'checksum': 'sha1$069b23dada099e20700c6b10471750e3cba6ee86',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.seg'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.seg'),
                                       'size': 647},
                                      {'basename': 'Sample1.qc.txt',
                                       'checksum': 'sha1$8f014a24b3e53c3801be6a6751de2dc089f79f69',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.qc.txt'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.qc.txt'),
                                       'size': 1224}],
                          'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2'),
                          'path': os.path.join(output_dir,'facets/Sample1.Sample2')
                            }]
                },
                'facets_failed_pairs': [],
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
                            'checksum': 'sha1$7d2bb282e74ff6a5d41b334ded689f9336722702',
                            'size': 132,
                            'path': os.path.join(output_dir, 'portal/meta_clinical_sample.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_patient.txt'),
                            'basename': 'data_clinical_patient.txt',
                            'class': 'File',
                            'checksum': 'sha1$cac1377a45dfc316266697a21df87801883127b5',
                            'size': 91,
                            'path': os.path.join(output_dir, 'portal/data_clinical_patient.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_sample.txt'),
                            'basename': 'data_clinical_sample.txt',
                            'class': 'File',
                            'checksum': 'sha1$678e003ed09df4f63a9a3a013623b7a18b6c9359',
                            'size': 1401,
                            'path': os.path.join(output_dir, 'portal/data_clinical_sample.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/meta_study.txt'),
                            'basename': 'meta_study.txt',
                            'class': 'File',
                            'checksum': 'sha1$998a850a03828cf4c235583dced9751ba75c9ab1',
                            'size': 128,
                            'path': os.path.join(output_dir, 'portal/meta_study.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_patient.txt'),
                            'basename': 'meta_clinical_patient.txt',
                            'class': 'File',
                            'checksum': 'sha1$e1f0b7786dd10af608df5178ff4b1a0b7a191a38',
                            'size': 134,
                            'path': os.path.join(output_dir, 'portal/meta_clinical_patient.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/meta_CNA.txt'),
                            'basename': 'meta_CNA.txt',
                            'class': 'File',
                            'checksum': 'sha1$93367ae36cae5e1a53b25e5bb02731e8b113251b',
                            'size': 262,
                            'path': os.path.join(output_dir, 'portal/meta_CNA.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/meta_fusions.txt'),
                            'basename': 'meta_fusions.txt',
                            'class': 'File',
                            'checksum': 'sha1$f23883a435cd6bcebb05b06ef6004ba3be413828',
                            'size': 219,
                            'path': os.path.join(output_dir, 'portal/meta_fusions.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/meta_mutations_extended.txt'),
                            'basename': 'meta_mutations_extended.txt',
                            'class': 'File',
                            'checksum': 'sha1$3be01b7a638dc3bb4f229038cbecd9b3cafc7410',
                            'size': 245,
                            'path': os.path.join(output_dir, 'portal/meta_mutations_extended.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt'),
                            'basename': 'Proj_08390_G_meta_cna_hg19_seg.txt',
                            'class': 'File',
                            'checksum': 'sha1$fd39b2b94b555902fb9b313a5d06cb7ad3a305eb',
                            'size': 192,
                            'path': os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.txt'),
                            'basename': 'data_CNA.txt',
                            'class': 'File',
                            'checksum': 'sha1$df07ff6b94392c470e2c8c47ca6e35cc2f6cc791',
                            'size': 2343,
                            'path': os.path.join(output_dir, 'portal/data_CNA.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                            'basename': 'data_CNA.ascna.txt',
                            'class': 'File',
                            'checksum': 'sha1$64710614a20b754dbbf863369ed18d05747166a3',
                            'size': 2759,
                            'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                            'basename': 'data_mutations_extended.txt',
                            'class': 'File',
                            'checksum': 'sha1$f4575aaf80b6b4e22cbfed5e7743d3285f7a9a88',
                            'size': 5623,
                            'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                            'basename': 'Proj_08390_G_data_cna_hg19.seg',
                            'class': 'File',
                            'checksum': 'sha1$e99e98d33598c423448cba44a37da712c9c933c5',
                            'size': 525,
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
                                'checksum': 'sha1$3d3d54b5cbd98e751adf926d07884eeaed94aaef',
                                'size': 188,
                                'path': os.path.join(output_dir, 'portal/case_lists/cases_all.txt')},
                                {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt'),
                                'basename': 'cases_cnaseq.txt',
                                'class': 'File',
                                'checksum': 'sha1$2a43c8ed31f69c12cf9a14044bab1aaa1b6ec7ee',
                                'size': 268,
                                'path': os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt')},
                                {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cna.txt'),
                                'basename': 'cases_cna.txt',
                                'class': 'File',
                                'checksum': 'sha1$117563761d7cb4690dea363f40b77ef1a0444fc1',
                                'size': 200,
                                'path': os.path.join(output_dir, 'portal/case_lists/cases_cna.txt')},
                                {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt'),
                                'basename': 'cases_sequenced.txt',
                                'class': 'File',
                                'checksum': 'sha1$78d1d0acb6928349b4448a73f2cdf1ed4edcb1a3',
                                'size': 213,
                                'path': os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt')}
                            ]
                        }
                    ]
                }
            }
        self.maxDiff = None
        self.assertDictEqual(output_json, expected_output)
        comments, mutations = self.load_mutations(os.path.join(output_dir, 'analysis', 'demo.muts.maf'))
        self.assertEqual(len(mutations), 22)
        comments, mutations = self.load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 17)

        # load the data_CNA.txt file
        path = os.path.join(output_dir, 'portal/data_CNA.txt') # renamed from the data_CNA.scna.txt file ...
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample1']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample1']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'portal/data_clinical_sample.txt')
        table_reader = TableReader(path)
        comments = table_reader.comment_lines
        fieldnames = table_reader.get_fieldnames()
        records = [ rec for rec in table_reader.read() ]

        expected_comments = [
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tgenome_doubled\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tASCN_WGD\tCMO_TMB_SCORE\n', '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tgenome_doubled\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tASCN_WGD\tCMO_TMB_SCORE\n', '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\tSTRING\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\tSTRING\tNUMBER\n', '#1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\n'
        ]
        self.assertEqual(comments, expected_comments)

        tmbs = {}
        for record in records:
            tmbs[record['SAMPLE_ID']] = record['CMO_TMB_SCORE']

        expected_tmbs = {'Sample1': '0.0', 'Sample4': 'NA'}
        self.assertEqual(tmbs, expected_tmbs)


    def test_run_worflow_one_maf(self):
        """
        Test that the workflow works correctly when run with a single maf
        """
        data_clinical_file = os.path.join(self.DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(self.DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        self.input = {
            "assay_coverage": "10000000", # TODO: get this from an assay reference key
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
                "path": self.IMPACT_FILE,
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
                "path": self.DATA_SETS['Proj_08390_G']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": self.KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "mutation_svs_maf_files": [
                {
                    "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.maf"),
                    "class": "File"
                }
            ],
            "pairs": [
                {
                    "pair_maf": {
                        "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                        "class": "File"
                    },
                    "snp_pileup": {
                         "path": os.path.join(self.DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads.dat.gz"),
                         "class": "File"
                    },
                    "pair_id": "Sample1.Sample2",
                    "tumor_id": "Sample1",
                    "normal_id": "Sample2"
                }
            ]
        }

        output_json, output_dir = self.run_cwl()

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
                        'checksum': 'sha1$890129aa22c2e6b0e6568d4799c628748fdc5139',
                        'size': 85269,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')
                    },
                    {
                        'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                        'basename': 'Proj_08390_G.muts.maf',
                        'class': 'File',
                        'checksum': 'sha1$5bea5122ba7ff08c40a41364decdf4be6dcff2c4',
                        'size': 39381,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')
                    },
                    {
                        'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf'),
                        'basename': 'Proj_08390_G.muts.share.maf',
                        'class': 'File',
                        'checksum': 'sha1$7716a2274203abe8039c99d0722b92a675ed38b7',
                        'size': 7288,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf')
                    },
                    {
                        'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt'),
                        'basename': 'Proj_08390_G.seg.cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$76e2c8a0a9c4200e26b8fc4a60701bd7b1b86560',
                        'size': 1459,
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
                'facets_dir': {
                    'location': 'file://' + os.path.join(output_dir, 'facets'),
                    'basename': 'facets',
                    'class': 'Directory',
                    'path': os.path.join(output_dir, 'facets'),
                    'listing': [{'basename': 'Sample1.Sample2',
                          'class': 'Directory',
                          'listing': [{'basename': 'Sample1.Sample2_hisens.ccf.portal.maf',
                                       'checksum': 'sha1$fa311e483303fbda7778622af7be57c7d62aac76',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.Sample2_hisens.ccf.portal.maf'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.Sample2_hisens.ccf.portal.maf'),
                                       'size': 19822273},
                                      {'basename': 'Sample1.arm_level.txt',
                                       'checksum': 'sha1$4f1d6620f3c613e40a147338739d3267e891dda7',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.arm_level.txt'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.arm_level.txt'),
                                       'size': 1986},
                                      {'basename': 'Sample1.txt',
                                       'checksum': 'sha1$8abe9bed5d1cb807b451c0647f658394d3f96708',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.txt'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.txt'),
                                       'size': 560},
                                      {'basename': 'Sample1.gene_level.txt',
                                       'checksum': 'sha1$6b31bff7758586ec39c94a0c2ddb8201eda5e0ce',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.gene_level.txt'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.gene_level.txt'),
                                       'size': 252793},
                                      {'basename': 'Sample1_hisens.cncf.txt',
                                       'checksum': 'sha1$a183b43abb9b7bb8f67ec6a0411539643d1275fe',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.cncf.txt'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.cncf.txt'),
                                       'size': 5774},
                                      {'basename': 'Sample1_hisens.rds',
                                       'checksum': 'sha1$d46031ce35efd7341146cac0db809c8cb75c97c2',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.rds'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.rds'),
                                       'size': 550930},
                                      {'basename': 'Sample1_hisens.seg',
                                       'checksum': 'sha1$1e0f1ad1cc751d53f02ad19717921027c3c96dfb',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.seg'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.seg'),
                                       'size': 1786},
                                      {'basename': 'Sample1_purity.rds',
                                       'checksum': 'sha1$265c7b94df391ae12e3bc1cc3ec4e763629a4f0a',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.rds'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.rds'),
                                       'size': 550388},
                                      {'basename': 'Sample1_purity.seg',
                                       'checksum': 'sha1$36904773ac5f194eac7f816a5fa9b055931bb48c',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.seg'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.seg'),
                                       'size': 1284},
                                      {'basename': 'Sample1.qc.txt',
                                       'checksum': 'sha1$6fdefe96a23c0bc52050c46a16625fca603e4eb6',
                                       'class': 'File',
                                       'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.qc.txt'),
                                       'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.qc.txt'),
                                       'size': 1340}],
                          'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2'),
                          'path': os.path.join(output_dir,'facets/Sample1.Sample2')
                            }]
                },
                'facets_failed_pairs': [],
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
                            'checksum': 'sha1$1d89d39e40c65f607dd71eedb21f93d935a44409',
                            'size': 8758,
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
                            'checksum': 'sha1$dca846b4f25ca20bca814fa15351de8d322d3439',
                            'size': 5319,
                            'path': os.path.join(output_dir, 'portal/data_CNA.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                            'basename': 'data_CNA.ascna.txt',
                            'class': 'File',
                            'checksum': 'sha1$a35bd02511775696127fb6bfcbd70ee75e51faba',
                            'size': 6167,
                            'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                            'basename': 'data_mutations_extended.txt',
                            'class': 'File',
                            'checksum': 'sha1$f98179c9fc244b91bd2521cc2dca7e611d5a8c9a',
                            'size': 5657,
                            'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                            'basename': 'Proj_08390_G_data_cna_hg19.seg',
                            'class': 'File',
                            'checksum': 'sha1$76e2c8a0a9c4200e26b8fc4a60701bd7b1b86560',
                            'size': 1459,
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
        comments, mutations = self.load_mutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'))
        self.assertEqual(len(mutations), 22)
        comments, mutations = self.load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 17)

        # load the data_CNA.txt file
        path = os.path.join(output_dir, 'portal/data_CNA.txt') # renamed from the data_CNA.scna.txt file ...
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample1']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample1']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'portal/data_clinical_sample.txt')
        table_reader = TableReader(path)
        comments = table_reader.comment_lines
        fieldnames = table_reader.get_fieldnames()
        records = [ rec for rec in table_reader.read() ]

        expected_comments = [
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tgenome_doubled\tASCN_WGD\tCMO_TMB_SCORE\n',
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tgenome_doubled\tASCN_WGD\tCMO_TMB_SCORE\n',
        '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\tSTRING\tSTRING\tNUMBER\n',
        '#1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\n'
        ]
        self.assertEqual(comments, expected_comments)

        tmbs = {}
        for record in records:
            tmbs[record['SAMPLE_ID']] = record['CMO_TMB_SCORE']

        expected_tmbs = {
        'Sample46': 'NA', 'Sample44': 'NA', 'Sample80': 'NA', 'Sample20': 'NA', 'Sample38': 'NA', 'Sample26': 'NA', 'Sample94': 'NA', 'Sample48': 'NA', 'Sample68': 'NA', 'Sample90': 'NA', 'Sample18': 'NA', 'Sample54': 'NA', 'Sample52': 'NA', 'Sample86': 'NA', 'Sample30': 'NA', 'Sample78': 'NA', 'Sample84': 'NA', 'Sample82': 'NA', 'Sample6': 'NA', 'Sample96': 'NA', 'Sample72': 'NA', 'Sample56': 'NA', 'Sample64': 'NA', 'Sample58': 'NA', 'Sample92': 'NA', 'Sample62': 'NA', 'Sample8': 'NA', 'Sample24': 'NA', 'Sample12': 'NA', 'Sample16': 'NA', 'Sample88': 'NA', 'Sample22': 'NA', 'Sample42': 'NA', 'Sample76': 'NA', 'Sample28': 'NA', 'Sample74': 'NA', 'Sample50': 'NA', 'Sample60': 'NA', 'Sample10': 'NA', 'Sample36': 'NA', 'Sample34': 'NA', 'Sample40': 'NA', 'Sample66': 'NA', 'Sample14': 'NA', 'Sample32': 'NA', 'Sample70': 'NA', 'Sample4': 'NA', 'Sample1': '0.0000000000013'
        }
        self.assertEqual(tmbs, expected_tmbs)

    def test_run_worflow_two_mafs(self):
        """
        Test that the workflow works correctly when run with two maf files
        """
        data_clinical_file = os.path.join(self.DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(self.DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        self.input = {
            "assay_coverage": "10000000", # TODO: get this from an assay reference key
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
                "path": self.IMPACT_FILE,
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
                "path": self.DATA_SETS['Proj_08390_G']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": self.KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "mutation_svs_maf_files": [
                {
                    "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.maf"),
                    "class": "File"
                },
                {
                    "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.svs.pass.vep.maf"),
                    "class": "File"
                }
            ],
            "pairs": [
                {
                    "pair_maf": {
                        "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                        "class": "File"
                    },
                    "snp_pileup": {
                         "path": os.path.join(self.DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads.dat.gz"),
                         "class": "File"
                    },
                    "pair_id": "Sample1.Sample2",
                    "tumor_id": "Sample1",
                    "normal_id": "Sample2"
                },
                {
                    "pair_maf": {
                        "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.muts.maf"),
                        "class": "File"
                    },
                    "snp_pileup": {
                         "path": os.path.join(self.DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample3.rg.md.abra.printreads__Sample4.rg.md.abra.printreads.dat.gz"),
                         "class": "File"
                    },
                    "pair_id": "Sample4.Sample3",
                    "tumor_id": "Sample4",
                    "normal_id": "Sample3"
                }
            ]
        }

        output_json, output_dir = self.run_cwl()

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
                    'checksum': 'sha1$adf367afca81f3bb6e75272c86c38c4ad20246a9',
                    'size': 170207,
                    'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')},
                    {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                    'basename': 'Proj_08390_G.muts.maf',
                    'class': 'File',
                    'checksum': 'sha1$213939883b558f956cbe007909c6d68e60a3add2',
                    'size': 63615,
                    'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')},
                    {
                        'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf'),
                        'basename': 'Proj_08390_G.muts.share.maf',
                        'class': 'File',
                        'checksum': 'sha1$2dafe49600cac4fbc00e230af45ee2eb258aa6ac',
                        'size': 10795,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf')
                    },
                    {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt'),
                    'basename': 'Proj_08390_G.seg.cna.txt',
                    'class': 'File',
                    'checksum': 'sha1$8995133edb5b4547371b0fdea9925848eb480fbc',
                    'size': 2703,
                    'path': os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt')},
                    {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf'),
                    'basename': 'Proj_08390_G.svs.maf',
                    'class': 'File',
                    'checksum': 'sha1$5c2a63fc01980550108e58079a8b689d53c97d8c',
                    'size': 35595,
                    'path': os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf')}
                ]
            },
            'facets_dir': {'basename': 'facets',
              'class': 'Directory',
              'listing': [{'basename': 'Sample1.Sample2',
                           'class': 'Directory',
                           'listing': [{'basename': 'Sample1.Sample2_hisens.ccf.portal.maf',
                                        'checksum': 'sha1$fa311e483303fbda7778622af7be57c7d62aac76',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.Sample2_hisens.ccf.portal.maf'),
                                        'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.Sample2_hisens.ccf.portal.maf'),
                                        'size': 19822273},
                                       {'basename': 'Sample1.arm_level.txt',
                                        'checksum': 'sha1$4f1d6620f3c613e40a147338739d3267e891dda7',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.arm_level.txt'),
                                        'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.arm_level.txt'),
                                        'size': 1986},
                                       {'basename': 'Sample1.txt',
                                        'checksum': 'sha1$8abe9bed5d1cb807b451c0647f658394d3f96708',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.txt'),
                                        'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.txt'),
                                        'size': 560},
                                       {'basename': 'Sample1.gene_level.txt',
                                        'checksum': 'sha1$6b31bff7758586ec39c94a0c2ddb8201eda5e0ce',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.gene_level.txt'),
                                        'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.gene_level.txt'),
                                        'size': 252793},
                                       {'basename': 'Sample1_hisens.cncf.txt',
                                        'checksum': 'sha1$a183b43abb9b7bb8f67ec6a0411539643d1275fe',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.cncf.txt'),
                                        'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.cncf.txt'),
                                        'size': 5774},
                                       {'basename': 'Sample1_hisens.rds',
                                        'checksum': 'sha1$d46031ce35efd7341146cac0db809c8cb75c97c2',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.rds'),
                                        'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.rds'),
                                        'size': 550930},
                                       {'basename': 'Sample1_hisens.seg',
                                        'checksum': 'sha1$1e0f1ad1cc751d53f02ad19717921027c3c96dfb',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.seg'),
                                        'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.seg'),
                                        'size': 1786},
                                       {'basename': 'Sample1_purity.rds',
                                        'checksum': 'sha1$265c7b94df391ae12e3bc1cc3ec4e763629a4f0a',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.rds'),
                                        'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.rds'),
                                        'size': 550388},
                                       {'basename': 'Sample1_purity.seg',
                                        'checksum': 'sha1$36904773ac5f194eac7f816a5fa9b055931bb48c',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.seg'),
                                        'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.seg'),
                                        'size': 1284},
                                       {'basename': 'Sample1.qc.txt',
                                        'checksum': 'sha1$6fdefe96a23c0bc52050c46a16625fca603e4eb6',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.qc.txt'),
                                        'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.qc.txt'),
                                        'size': 1340}],
                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2'),
                           'path': os.path.join(output_dir,'facets/Sample1.Sample2')},
                          {'basename': 'Sample4.Sample3',
                           'class': 'Directory',
                           'listing': [{'basename': 'Sample4.Sample3_hisens.ccf.portal.maf',
                                        'checksum': 'sha1$86af92550b3e4252f950223f997579bfe651b65f',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.Sample3_hisens.ccf.portal.maf'),
                                        'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.Sample3_hisens.ccf.portal.maf'),
                                        'size': 1136054},
                                       {'basename': 'Sample4.arm_level.txt',
                                        'checksum': 'sha1$09b97c42a9875d642dae4d030712264f1b92b550',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.arm_level.txt'),
                                        'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.arm_level.txt'),
                                        'size': 1424},
                                       {'basename': 'Sample4.txt',
                                        'checksum': 'sha1$1bb6cf1ea21bc774037913633e5e06772d791b35',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.txt'),
                                        'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.txt'),
                                        'size': 602},
                                       {'basename': 'Sample4.gene_level.txt',
                                        'checksum': 'sha1$4a1f580aebf86a3a2ce4d5d27bc1caa3655e8788',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.gene_level.txt'),
                                        'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.gene_level.txt'),
                                        'size': 144319},
                                       {'basename': 'Sample4_hisens.cncf.txt',
                                        'checksum': 'sha1$a51eac1442aac591412b79305f003d464322140e',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.cncf.txt'),
                                        'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.cncf.txt'),
                                        'size': 4802},
                                       {'basename': 'Sample4_hisens.rds',
                                        'checksum': 'sha1$33286048372616b019fcf52587b25bbbf5d1ae08',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.rds'),
                                        'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.rds'),
                                        'size': 186732},
                                       {'basename': 'Sample4_hisens.seg',
                                        'checksum': 'sha1$ea32de058625243fdc9f48d9cbab746444549e2e',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.seg'),
                                        'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.seg'),
                                        'size': 1584},
                                       {'basename': 'Sample4_purity.rds',
                                        'checksum': 'sha1$1eda8dc016cf8322cab1190ea4f723ecb33c4440',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_purity.rds'),
                                        'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_purity.rds'),
                                        'size': 186513},
                                       {'basename': 'Sample4_purity.seg',
                                        'checksum': 'sha1$6be6775133173d609e176f16a70c5a3c2aeb9027',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_purity.seg'),
                                        'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_purity.seg'),
                                        'size': 1276},
                                       {'basename': 'Sample4.qc.txt',
                                        'checksum': 'sha1$0c17bd63f43a6726682dec93a58cbcb3ccfd20e1',
                                        'class': 'File',
                                        'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.qc.txt'),
                                        'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.qc.txt'),
                                        'size': 1330}],
                           'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3'),
                           'path': os.path.join(output_dir,'facets/Sample4.Sample3')}],
              'location': 'file://' + os.path.join(output_dir,'facets'),
              'path': os.path.join(output_dir,'facets')
            },
            'facets_failed_pairs': [],
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
                    'checksum': 'sha1$013055a53721fece2488a88c16751f8a8dd26901',
                    'size': 8769,
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
                    'checksum': 'sha1$4795ba3d6c30a028d21bafbb55227ede294b61d5',
                    'size': 6727,
                    'path': os.path.join(output_dir, 'portal/data_CNA.txt')},
                    {'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                    'basename': 'data_CNA.ascna.txt',
                    'class': 'File',
                    'checksum': 'sha1$1081a3ab9830e0d1eaddafb45e93956264fdfba9',
                    'size': 8796,
                    'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')},
                    {'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                    'basename': 'data_mutations_extended.txt',
                    'class': 'File',
                    'checksum': 'sha1$e688086833d015e38e4dea2e207552190bae5d21',
                    'size': 8382,
                    'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')},
                    {'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                    'basename': 'Proj_08390_G_data_cna_hg19.seg',
                    'class': 'File',
                    'checksum': 'sha1$8995133edb5b4547371b0fdea9925848eb480fbc',
                    'size': 2703,
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
        comments, mutations = self.load_mutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'))
        self.assertEqual(len(mutations), 34)
        comments, mutations = self.load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 27)

        # the clonality column needs to have been added in the workflow output
        for mut in mutations:
            self.assertTrue('ASCN.CLONAL' in mut)

        # load the data_CNA.txt file
        path = os.path.join(output_dir, 'portal/data_CNA.txt') # renamed from the data_CNA.scna.txt file ...
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample1', 'Sample4']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample1', 'Sample4']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'portal/data_clinical_sample.txt')
        table_reader = TableReader(path)
        comments = table_reader.comment_lines
        fieldnames = table_reader.get_fieldnames()
        records = [ rec for rec in table_reader.read() ]

        expected_comments = [
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tgenome_doubled\tASCN_WGD\tCMO_TMB_SCORE\n',
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tgenome_doubled\tASCN_WGD\tCMO_TMB_SCORE\n',
        '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\tSTRING\tSTRING\tNUMBER\n',
        '#1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\n'
        ]
        self.assertEqual(comments, expected_comments)

        tmbs = {}
        for record in records:
            tmbs[record['SAMPLE_ID']] = record['CMO_TMB_SCORE']

        expected_tmbs = {
        'Sample46': 'NA', 'Sample44': 'NA', 'Sample80': 'NA', 'Sample20': 'NA', 'Sample38': 'NA', 'Sample26': 'NA', 'Sample94': 'NA', 'Sample48': 'NA', 'Sample68': 'NA', 'Sample90': 'NA', 'Sample18': 'NA', 'Sample54': 'NA', 'Sample52': 'NA', 'Sample86': 'NA', 'Sample30': 'NA', 'Sample78': 'NA', 'Sample84': 'NA', 'Sample82': 'NA', 'Sample6': 'NA', 'Sample96': 'NA', 'Sample72': 'NA', 'Sample56': 'NA', 'Sample64': 'NA', 'Sample58': 'NA', 'Sample92': 'NA', 'Sample62': 'NA', 'Sample8': 'NA', 'Sample24': 'NA', 'Sample12': 'NA', 'Sample16': 'NA', 'Sample88': 'NA', 'Sample22': 'NA', 'Sample42': 'NA', 'Sample76': 'NA', 'Sample28': 'NA', 'Sample74': 'NA', 'Sample50': 'NA', 'Sample60': 'NA', 'Sample10': 'NA', 'Sample36': 'NA', 'Sample34': 'NA', 'Sample40': 'NA', 'Sample66': 'NA', 'Sample14': 'NA', 'Sample32': 'NA', 'Sample70': 'NA', 'Sample4': '0.0', 'Sample1': '0.0000000000013'
        }
        self.assertEqual(tmbs, expected_tmbs)

if __name__ == "__main__":
    unittest.main()
