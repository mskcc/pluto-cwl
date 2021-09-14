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
from pluto.settings import ENABLE_LARGE_TESTS, MICROSATELLITES_LIST
from pluto.serializer import OFile, ODir
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
        normal_bam = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample2.bam")
        tumor_bam = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample1.bam")
        microsatellites_file = self.DATA_SETS['demo']['microsatellites_file']

        self.input = {
            "assay_coverage": "1000", # TODO: get this from an assay reference key
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
            "microsatellites_file": {
                "path": microsatellites_file,
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
            ],
            # these must be in the same order as pairs
            "normal_bam_files": [
                {'class': 'File', 'path': normal_bam}
            ],
            "tumor_bam_files": [
                {'class': 'File', 'path': tumor_bam}
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_dir': ODir(name = 'analysis', dir = output_dir, items = [
                OFile(name = 'demo.gene.cna.txt', size = 26505, hash = 'a3527d2bf248c3121f093f7edfdf2605b12dfe71'),
                OFile(name = 'demo.muts.maf', size = 36333, hash = 'bbb1320d792c5a37f0d2afae82fdd3e67b2db182'),
                OFile(name = 'demo.muts.share.maf', size = 7288, hash = '7716a2274203abe8039c99d0722b92a675ed38b7'),
                OFile(name = 'demo.seg.cna.txt', size = 525, hash = 'e99e98d33598c423448cba44a37da712c9c933c5'),
                OFile(name = 'demo.svs.maf', size = 23603, hash = 'df420706bb5b772a79317843c0a01a3c88a9571d')
            ]),
                'facets_dir': ODir(name = 'facets', dir = output_dir, items = [
                    ODir(name = 'Sample1.Sample2', items = [
                        OFile(name = 'Sample1.Sample2_hisens.ccf.portal.maf', size = 62216, hash = 'e86566d281bc808e57ee912d02bdbda8179d7c0a'),
                        OFile(name = 'Sample1.arm_level.txt', size = 370, hash = '98704689aedd071919ca0b2e387462948cb803b4'),
                        OFile(name = 'Sample1.txt', size = 506, hash = 'ccdaded463980977bc4c90e539043c6efd1aef33'),
                        OFile(name = 'Sample1.gene_level.txt', size = 94824, hash = '75ad4054ccec1f7a47d1c7738730816191d713c3'),
                        OFile(name = 'Sample1_hisens.cncf.txt', size = 1385, hash = '738974eaa3e4278f40cdbbf487b87f8b57ee2ecb'),
                        OFile(name = 'Sample1_hisens.rds', size = 191089, hash = '13846db5dc75ca6b7d64cb69ff04fa64da9a88e7'),
                        OFile(name = 'Sample1_hisens.seg', size = 647, hash = '069b23dada099e20700c6b10471750e3cba6ee86'),
                        OFile(name = 'Sample1_purity.rds', size = 191195, hash = 'ded479ac6278f1cd7bf42bf186440d8a89415a02'),
                        OFile(name = 'Sample1_purity.seg', size = 647, hash = '069b23dada099e20700c6b10471750e3cba6ee86'),
                        OFile(name = 'Sample1.qc.txt', size = 1224, hash = '8f014a24b3e53c3801be6a6751de2dc089f79f69'),
                    ])
                ]),
                'facets_failed_pairs': [],
                'portal_dir': ODir(name = 'portal', dir = output_dir, items = [
                    OFile(name = 'meta_clinical_sample.txt', size = 132, hash = '7d2bb282e74ff6a5d41b334ded689f9336722702'),
                    OFile(name = 'data_clinical_patient.txt', size = 91, hash = 'cac1377a45dfc316266697a21df87801883127b5'),
                    OFile(name = 'data_clinical_sample.txt', size = 1523, hash = '40caf8f2180df9f8b7e9af8b6913f27274d3aa87'),
                    OFile(name = 'meta_study.txt', size = 128, hash = '998a850a03828cf4c235583dced9751ba75c9ab1'),
                    OFile(name = 'meta_clinical_patient.txt', size = 134, hash = 'e1f0b7786dd10af608df5178ff4b1a0b7a191a38'),
                    OFile(name = 'meta_CNA.txt', size = 262, hash = '93367ae36cae5e1a53b25e5bb02731e8b113251b'),
                    OFile(name = 'meta_fusions.txt', size = 219, hash = 'f23883a435cd6bcebb05b06ef6004ba3be413828'),
                    OFile(name = 'meta_mutations_extended.txt', size = 245, hash = '3be01b7a638dc3bb4f229038cbecd9b3cafc7410'),
                    OFile(name = 'Proj_08390_G_meta_cna_hg19_seg.txt', size = 192, hash = 'fd39b2b94b555902fb9b313a5d06cb7ad3a305eb'),
                    OFile(name = 'data_CNA.txt', size = 2343, hash = 'df07ff6b94392c470e2c8c47ca6e35cc2f6cc791'),
                    OFile(name = 'data_CNA.ascna.txt', size = 2759, hash = '64710614a20b754dbbf863369ed18d05747166a3'),
                    OFile(name = 'data_mutations_extended.txt', size = 5623, hash = 'f4575aaf80b6b4e22cbfed5e7743d3285f7a9a88'),
                    OFile(name = 'Proj_08390_G_data_cna_hg19.seg', size = 525, hash = 'e99e98d33598c423448cba44a37da712c9c933c5'),
                    OFile(name = 'data_fusions.txt', size = 99, hash = 'c16f763b248813fcdde76f7486f1ddc4e9856038'),
                    ODir(name = 'case_lists', items = [
                        OFile(name = 'cases_all.txt', size = 188, hash = '3d3d54b5cbd98e751adf926d07884eeaed94aaef'),
                        OFile(name = 'cases_cnaseq.txt', size = 268, hash = '2a43c8ed31f69c12cf9a14044bab1aaa1b6ec7ee'),
                        OFile(name = 'cases_cna.txt', size = 200, hash = '117563761d7cb4690dea363f40b77ef1a0444fc1'),
                        OFile(name = 'cases_sequenced.txt', size = 213, hash = '78d1d0acb6928349b4448a73f2cdf1ed4edcb1a3'),
                    ]),
                    OFile(name = 'report.html')
                ])
            }
        for item in output_json['portal_dir']['listing']:
            if item['basename'] == 'report.html':
                item.pop('checksum')
                item.pop('size')
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
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tgenome_doubled\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tASCN_WGD\tCMO_TMB_SCORE\tCMO_MSI_SCORE\tCMO_MSI_STATUS\n', '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tgenome_doubled\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tASCN_WGD\tCMO_TMB_SCORE\tCMO_MSI_SCORE\tCMO_MSI_STATUS\n', '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\tSTRING\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\n', '#1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t0\t0\n'
        ]
        self.assertEqual(comments, expected_comments)


        metrics = {}
        metrics['tmb'] = {}
        metrics['msi'] = {}
        metrics['msi_status'] = {}
        for record in records:
            sample_id = record['SAMPLE_ID']
            metrics['tmb'][sample_id] = record['CMO_TMB_SCORE']
            metrics['msi'][sample_id] = record['MSI_SCORE']
            metrics['msi_status'][sample_id] = record['MSI_STATUS']

        expected_tmbs = {'Sample1': '17000.0', 'Sample4': 'NA'}
        self.assertEqual(metrics['tmb'], expected_tmbs)

        expected_msis = {'Sample1': '28.00', 'Sample4': 'NA'}
        self.assertEqual(metrics['msi'], expected_msis)

        expected_msi_statuses = {'Sample1': 'Instable', 'Sample4': 'NA'}
        self.assertEqual(metrics['msi_status'], expected_msi_statuses)

    #################################################################
    #
    #  ########  ######## ##     ##  #######   #######
    #  ##     ## ##       ###   ### ##     ## ##     ##
    #  ##     ## ##       #### #### ##     ##        ##
    #  ##     ## ######   ## ### ## ##     ##  #######
    #  ##     ## ##       ##     ## ##     ## ##
    #  ##     ## ##       ##     ## ##     ## ##
    #  ########  ######## ##     ##  #######  #########
    #
    # https://patorjk.com/software/taag/#p=display&c=bash&f=Banner3&t=DEMO2
    #################################################################
    def test_demo_dataset2(self):
        """
        Test case for two small demo datasets
        """
        data_clinical_file = os.path.join(self.DATA_SETS['demo']['INPUTS_DIR'], "demo_sample_data_clinical.txt")
        sample_summary_file = os.path.join(self.DATA_SETS['demo']['QC_DIR'], "demo_SampleSummary.txt")
        microsatellites_file = self.DATA_SETS['demo']['microsatellites_file']

        mutation_svs_txt_file1 = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt")
        mutation_svs_maf1 = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.maf")
        pair_maf1 = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        snp_pileup1 = os.path.join(self.DATA_SETS['demo']['SNP_PILEUP_DIR'], "Sample1.Sample2.snp_pileup.gz")
        normal_bam1 = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample2.bam")
        tumor_bam1 = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample1.bam")

        mutation_svs_txt_file2 = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample4.Sample3.svs.pass.vep.portal.txt")
        mutation_svs_maf2 = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample4.Sample3.svs.pass.vep.maf")
        pair_maf2 = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample4.Sample3.muts.maf")
        snp_pileup2 = os.path.join(self.DATA_SETS['demo']['SNP_PILEUP_DIR'], "Sample4.Sample3.snp_pileup.gz")
        normal_bam2 = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample3.bam")
        tumor_bam2 = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample4.bam")

        self.input = {
            "assay_coverage": "1000", # TODO: get this from an assay reference key
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
            "microsatellites_file": {
                "path": microsatellites_file,
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
            "mutation_svs_txt_files": [
                { "path": mutation_svs_txt_file1, "class": "File" },
                { "path": mutation_svs_txt_file2, "class": "File" }
                 ],
            "mutation_svs_maf_files": [
                { "path": mutation_svs_maf1, "class": "File" } ,
                { "path": mutation_svs_maf2, "class": "File" } ,
                ],
            "pairs": [
                {
                    "pair_maf": { "path": pair_maf1, "class": "File" },
                    "snp_pileup": { "path": snp_pileup1, "class": "File" },
                    "pair_id": "Sample1.Sample2",
                    "tumor_id": "Sample1",
                    "normal_id": "Sample2"
                },
                {
                    "pair_maf": { "path": pair_maf2, "class": "File" },
                    "snp_pileup": { "path": snp_pileup2, "class": "File" },
                    "pair_id": "Sample4.Sample3",
                    "tumor_id": "Sample4",
                    "normal_id": "Sample3"
                },
            ],
            # these must be in the same order as pairs
            "normal_bam_files": [
                {'class': 'File', 'path': normal_bam1},
                {'class': 'File', 'path': normal_bam2}
            ],
            "tumor_bam_files": [
                {'class': 'File', 'path': tumor_bam1},
                {'class': 'File', 'path': tumor_bam2}
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_dir': ODir(name = 'analysis', dir = output_dir, items = [
                OFile(name = 'demo.gene.cna.txt', size = 41818, hash = 'f68fbc5697ef2437c185e8a2ceef926b023de3dc'),
                OFile(name = 'demo.muts.maf', size = 59074, hash = 'a8eec9ff26a8d3fe47bc301d34b1a475899f0b48'),
                OFile(name = 'demo.muts.share.maf', size = 10795, hash = '22f2aa10ff2080458fc895ab93d2e9707d634b32'),
                OFile(name = 'demo.seg.cna.txt', size = 878, hash = 'e1c85484cd13a88cea322de84cd35143c86a8ef3'),
                OFile(name = 'demo.svs.maf', size = 35595, hash = '5c2a63fc01980550108e58079a8b689d53c97d8c')
            ]),
                'facets_dir': ODir(name = 'facets', dir = output_dir, items = [
                    ODir(name = 'Sample1.Sample2', items = [
                        OFile(name = 'Sample1.Sample2_hisens.ccf.portal.maf', size = 62216, hash = 'e86566d281bc808e57ee912d02bdbda8179d7c0a'),
                        OFile(name = 'Sample1.arm_level.txt', size = 370, hash = '98704689aedd071919ca0b2e387462948cb803b4'),
                        OFile(name = 'Sample1.txt', size = 506, hash = 'ccdaded463980977bc4c90e539043c6efd1aef33'),
                        OFile(name = 'Sample1.gene_level.txt', size = 94824, hash = '75ad4054ccec1f7a47d1c7738730816191d713c3'),
                        OFile(name = 'Sample1_hisens.cncf.txt', size = 1385, hash = '738974eaa3e4278f40cdbbf487b87f8b57ee2ecb'),
                        OFile(name = 'Sample1_hisens.rds', size = 191089, hash = '13846db5dc75ca6b7d64cb69ff04fa64da9a88e7'),
                        OFile(name = 'Sample1_hisens.seg', size = 647, hash = '069b23dada099e20700c6b10471750e3cba6ee86'),
                        OFile(name = 'Sample1_purity.rds', size = 191195, hash = 'ded479ac6278f1cd7bf42bf186440d8a89415a02'),
                        OFile(name = 'Sample1_purity.seg', size = 647, hash = '069b23dada099e20700c6b10471750e3cba6ee86'),
                        OFile(name = 'Sample1.qc.txt', size = 1224, hash = '8f014a24b3e53c3801be6a6751de2dc089f79f69'),
                    ]),
                    ODir(name = 'Sample4.Sample3', items = [
                        OFile(name = 'Sample4.Sample3_hisens.ccf.portal.maf', size = 24450, hash = 'dd7be7c4a6256ad6a1d6a168fe24e5b347a0bb31'),
                        OFile(name = 'Sample4.arm_level.txt', size = 231, hash = 'e062db543739dbd96466a99a3237de1e27c1fcc3'),
                        OFile(name = 'Sample4.txt', size = 508, hash = 'b20d6ec5d51ef3cd8564ae1f3281c76db8b938a0'),
                        OFile(name = 'Sample4.gene_level.txt', size = 35806, hash = '2cc5c01ef531dff22610daf156f1f31e669b59a0'),
                        OFile(name = 'Sample4_hisens.cncf.txt', size = 1056, hash = '946087e73a21c3cccf2f354e533b9a13bbb8e732'),
                        OFile(name = 'Sample4_hisens.rds', size = 55045, hash = 'e2d0bd96c6d63ee70422a6976f2e72d89f1c2e80'),
                        OFile(name = 'Sample4_hisens.seg', size = 488, hash = 'e6df130c57ca594578f9658e589cfafc8f40a56c'),
                        OFile(name = 'Sample4_purity.rds', size = 55150, hash = '5eea298ef9bfe52096dfd18825ad47adb52a1f99'),
                        OFile(name = 'Sample4_purity.seg', size = 488, hash = 'e6df130c57ca594578f9658e589cfafc8f40a56c'),
                        OFile(name = 'Sample4.qc.txt', size = 1230, hash = 'd862c00273178e537eeeccc69db15cfa20e811b2'),
                    ]),
                ]),
                'facets_failed_pairs': [],
                'portal_dir': ODir(name = 'portal', dir = output_dir, items = [
                    OFile(name = 'meta_clinical_sample.txt', size = 132, hash = '7d2bb282e74ff6a5d41b334ded689f9336722702'),
                    OFile(name = 'data_clinical_patient.txt', size = 91, hash = 'cac1377a45dfc316266697a21df87801883127b5'),
                    OFile(name = 'data_clinical_sample.txt', size = 1546, hash = 'a975d9be1ade1185a2d78ed3be0d45ffb584be52'),
                    OFile(name = 'meta_study.txt', size = 128, hash = '998a850a03828cf4c235583dced9751ba75c9ab1'),
                    OFile(name = 'meta_clinical_patient.txt', size = 134, hash = 'e1f0b7786dd10af608df5178ff4b1a0b7a191a38'),
                    OFile(name = 'meta_CNA.txt', size = 262, hash = '93367ae36cae5e1a53b25e5bb02731e8b113251b'),
                    OFile(name = 'meta_fusions.txt', size = 219, hash = 'f23883a435cd6bcebb05b06ef6004ba3be413828'),
                    OFile(name = 'meta_mutations_extended.txt', size = 245, hash = '3be01b7a638dc3bb4f229038cbecd9b3cafc7410'),
                    OFile(name = 'Proj_08390_G_meta_cna_hg19_seg.txt', size = 192, hash = 'fd39b2b94b555902fb9b313a5d06cb7ad3a305eb'),
                    OFile(name = 'data_CNA.txt', size = 3808, hash = 'ddc7acafdb61a3ea2aa1ec02e01df7bcc1691267'),
                    OFile(name = 'data_CNA.ascna.txt', size = 5104, hash = 'e1570ff3ae4202b99322914b5edb0579c520adaa'),
                    OFile(name = 'data_mutations_extended.txt', size = 8338, hash = '581fde39a445050f288a510beb4623eaf5346346'),
                    OFile(name = 'Proj_08390_G_data_cna_hg19.seg', size = 878, hash = 'e1c85484cd13a88cea322de84cd35143c86a8ef3'),
                    OFile(name = 'data_fusions.txt', size = 99, hash = 'c16f763b248813fcdde76f7486f1ddc4e9856038'),
                    ODir(name = 'case_lists', items = [
                        OFile(name = 'cases_all.txt', size = 188, hash = '3d3d54b5cbd98e751adf926d07884eeaed94aaef'),
                        OFile(name = 'cases_cnaseq.txt', size = 268, hash = '2a43c8ed31f69c12cf9a14044bab1aaa1b6ec7ee'),
                        OFile(name = 'cases_cna.txt', size = 200, hash = '117563761d7cb4690dea363f40b77ef1a0444fc1'),
                        OFile(name = 'cases_sequenced.txt', size = 213, hash = '78d1d0acb6928349b4448a73f2cdf1ed4edcb1a3'),
                    ]),
                    OFile(name = 'report.html')
                ])
            }

        for item in output_json['portal_dir']['listing']:
            if item['basename'] == 'report.html':
                item.pop('checksum')
                item.pop('size')

        self.maxDiff = None
        self.assertDictEqual(output_json, expected_output)
        comments, mutations = self.load_mutations(os.path.join(output_dir, 'analysis', 'demo.muts.maf'))
        self.assertEqual(len(mutations), 34)
        comments, mutations = self.load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 27)

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
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tgenome_doubled\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tASCN_WGD\tCMO_TMB_SCORE\tCMO_MSI_SCORE\tCMO_MSI_STATUS\n', '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tgenome_doubled\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tASCN_WGD\tCMO_TMB_SCORE\tCMO_MSI_SCORE\tCMO_MSI_STATUS\n', '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\tSTRING\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\n', '#1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t0\t0\n'
        ]
        self.assertEqual(comments, expected_comments)

        metrics = {}
        metrics['tmb'] = {}
        metrics['msi'] = {}
        metrics['msi_status'] = {}
        for record in records:
            sample_id = record['SAMPLE_ID']
            metrics['tmb'][sample_id] = record['CMO_TMB_SCORE']
            metrics['msi'][sample_id] = record['MSI_SCORE']
            metrics['msi_status'][sample_id] = record['MSI_STATUS']

        expected_tmbs = {'Sample1': '17000.0', 'Sample4': '9000.0'}
        self.assertEqual(metrics['tmb'], expected_tmbs)

        expected_msis = {'Sample1': '28.00', 'Sample4': '11.76'}
        self.assertEqual(metrics['msi'], expected_msis)

        expected_msi_statuses = {'Sample1': 'Instable', 'Sample4': 'Instable'}
        self.assertEqual(metrics['msi_status'], expected_msi_statuses)

    #################################################################
    #
    #  ######## ##     ## ##       ##             ##
    #  ##       ##     ## ##       ##           ####
    #  ##       ##     ## ##       ##             ##
    #  ######   ##     ## ##       ##             ##
    #  ##       ##     ## ##       ##             ##
    #  ##       ##     ## ##       ##             ##
    #  ##        #######  ######## ########     ######
    #
    #################################################################
    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_run_worflow_one_maf(self):
        """
        Test that the workflow works correctly when run with a single maf
        """
        data_clinical_file = os.path.join(self.DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(self.DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        microsatellites_file = self.DATA_SETS['demo']['microsatellites_file']
        tumor_bam = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample1.rg.md.abra.printreads.bam")
        normal_bam = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample2.rg.md.abra.printreads.bam")

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
            "microsatellites_file": {
                "path": microsatellites_file,
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
            ],
            # these must be in the same order as pairs
            "normal_bam_files": [
                {'class': 'File', 'path': normal_bam}
            ],
            "tumor_bam_files": [
                {'class': 'File', 'path': tumor_bam}
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_dir': ODir(name = 'analysis', dir = output_dir, items = [
                OFile(name = 'Proj_08390_G.gene.cna.txt', size = 85269, hash = '890129aa22c2e6b0e6568d4799c628748fdc5139'),
                OFile(name = 'Proj_08390_G.muts.maf', size = 39381, hash = '5bea5122ba7ff08c40a41364decdf4be6dcff2c4'),
                OFile(name = 'Proj_08390_G.muts.share.maf', size = 7288, hash = '7716a2274203abe8039c99d0722b92a675ed38b7'),
                OFile(name = 'Proj_08390_G.seg.cna.txt', size = 1459, hash = '76e2c8a0a9c4200e26b8fc4a60701bd7b1b86560'),
                OFile(name = 'Proj_08390_G.svs.maf', size = 23603, hash = 'df420706bb5b772a79317843c0a01a3c88a9571d')
            ]),
                'facets_dir': ODir(name = 'facets', dir = output_dir, items = [
                    ODir(name = 'Sample1.Sample2', items = [
                        OFile(name = 'Sample1.Sample2_hisens.ccf.portal.maf', size = 19822273, hash = 'fa311e483303fbda7778622af7be57c7d62aac76'),
                        OFile(name = 'Sample1.arm_level.txt', size = 1986, hash = '4f1d6620f3c613e40a147338739d3267e891dda7'),
                        OFile(name = 'Sample1.txt', size = 560, hash = '8abe9bed5d1cb807b451c0647f658394d3f96708'),
                        OFile(name = 'Sample1.gene_level.txt', size = 252793, hash = '6b31bff7758586ec39c94a0c2ddb8201eda5e0ce'),
                        OFile(name = 'Sample1_hisens.cncf.txt', size = 5774, hash = 'a183b43abb9b7bb8f67ec6a0411539643d1275fe'),
                        OFile(name = 'Sample1_hisens.rds', size = 550930, hash = 'd46031ce35efd7341146cac0db809c8cb75c97c2'),
                        OFile(name = 'Sample1_hisens.seg', size = 1786, hash = '1e0f1ad1cc751d53f02ad19717921027c3c96dfb'),
                        OFile(name = 'Sample1_purity.rds', size = 550388, hash = '265c7b94df391ae12e3bc1cc3ec4e763629a4f0a'),
                        OFile(name = 'Sample1_purity.seg', size = 1284, hash = '36904773ac5f194eac7f816a5fa9b055931bb48c'),
                        OFile(name = 'Sample1.qc.txt', size = 1340, hash = '6fdefe96a23c0bc52050c46a16625fca603e4eb6'),
                    ]),
                ]),
                'facets_failed_pairs': [],
                'portal_dir': ODir(name = 'portal', dir = output_dir, items = [
                    OFile(name = 'meta_clinical_sample.txt', size = 140, hash = '4c567d81c3b17a76c324fd3e2f73793a6e804f65'),
                    OFile(name = 'data_clinical_patient.txt', size = 643, hash = '9417dcabddd6ab2cbe98167bccd9b9e4fa182562'),
                    OFile(name = 'data_clinical_sample.txt', size = 9141, hash = 'b700cf1afea9e4223fd160d49681815f69c61401'),
                    OFile(name = 'meta_study.txt', size = 152, hash = '2b0a5fd1a97329adf7c3b1596c84cd6567059a95'),
                    OFile(name = 'meta_clinical_patient.txt', size = 142, hash = '9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f'),
                    OFile(name = 'meta_CNA.txt', size = 270, hash = 'a9bf16f6a0490b19e611e8814b85f7bf1d52417a'),
                    OFile(name = 'meta_fusions.txt', size = 227, hash = '77649e888bafc6a4ed61261d1c46d2f238e1c32b'),
                    OFile(name = 'meta_mutations_extended.txt', size = 253, hash = 'fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d'),
                    OFile(name = 'Proj_08390_G_meta_cna_hg19_seg.txt', size = 200, hash = '59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e'),
                    OFile(name = 'data_CNA.txt', size = 5319, hash = 'dca846b4f25ca20bca814fa15351de8d322d3439'),
                    OFile(name = 'data_CNA.ascna.txt', size = 6167, hash = 'a35bd02511775696127fb6bfcbd70ee75e51faba'),
                    OFile(name = 'data_mutations_extended.txt', size = 5657, hash = 'f98179c9fc244b91bd2521cc2dca7e611d5a8c9a'),
                    OFile(name = 'Proj_08390_G_data_cna_hg19.seg', size = 1459, hash = '76e2c8a0a9c4200e26b8fc4a60701bd7b1b86560'),
                    OFile(name = 'data_fusions.txt', size = 99, hash = 'c16f763b248813fcdde76f7486f1ddc4e9856038'),
                    ODir(name = 'case_lists', items = [
                        OFile(name = 'cases_all.txt', size = 616, hash = 'b9e43289cec5603b0886b5e8507c8d019387c125'),
                        OFile(name = 'cases_cnaseq.txt', size = 696, hash = 'b87e2da8dce0fddbadec348efe2986519b2a794b'),
                        OFile(name = 'cases_cna.txt', size = 628, hash = '053481a8299e9430117f8e45e081aa7ec21033a6'),
                        OFile(name = 'cases_sequenced.txt', size = 641, hash = 'ef9f5aef03c2527bf576470168660557ca1c7cc9'),
                    ]),
                    OFile(name = 'report.html')
                ])
            }

        for item in output_json['portal_dir']['listing']:
            if item['basename'] == 'report.html':
                item.pop('checksum')
                item.pop('size')

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
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tgenome_doubled\tASCN_WGD\tCMO_TMB_SCORE\tCMO_MSI_SCORE\tCMO_MSI_STATUS\n',
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tgenome_doubled\tASCN_WGD\tCMO_TMB_SCORE\tCMO_MSI_SCORE\tCMO_MSI_STATUS\n',
        '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\n',
        '#1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t0\t0\n'
        ]
        self.assertEqual(comments, expected_comments)

        metrics = {}
        metrics['tmb'] = {}
        metrics['msi'] = {}
        metrics['msi_status'] = {}
        for record in records:
            sample_id = record['SAMPLE_ID']
            metrics['tmb'][sample_id] = record['CMO_TMB_SCORE']
            metrics['msi'][sample_id] = record['MSI_SCORE']
            metrics['msi_status'][sample_id] = record['MSI_STATUS']

        expected_tmbs = {
        'Sample46': 'NA', 'Sample44': 'NA', 'Sample80': 'NA', 'Sample20': 'NA', 'Sample38': 'NA', 'Sample26': 'NA', 'Sample94': 'NA', 'Sample48': 'NA', 'Sample68': 'NA', 'Sample90': 'NA', 'Sample18': 'NA', 'Sample54': 'NA', 'Sample52': 'NA', 'Sample86': 'NA', 'Sample30': 'NA', 'Sample78': 'NA', 'Sample84': 'NA', 'Sample82': 'NA', 'Sample6': 'NA', 'Sample96': 'NA', 'Sample72': 'NA', 'Sample56': 'NA', 'Sample64': 'NA', 'Sample58': 'NA', 'Sample92': 'NA', 'Sample62': 'NA', 'Sample8': 'NA', 'Sample24': 'NA', 'Sample12': 'NA', 'Sample16': 'NA', 'Sample88': 'NA', 'Sample22': 'NA', 'Sample42': 'NA', 'Sample76': 'NA', 'Sample28': 'NA', 'Sample74': 'NA', 'Sample50': 'NA', 'Sample60': 'NA', 'Sample10': 'NA', 'Sample36': 'NA', 'Sample34': 'NA', 'Sample40': 'NA', 'Sample66': 'NA', 'Sample14': 'NA', 'Sample32': 'NA', 'Sample70': 'NA', 'Sample4': 'NA', 'Sample1': '47.5'
        }
        self.assertEqual(metrics['tmb'], expected_tmbs)

        expected_msis = {'Sample46': 'NA', 'Sample44': 'NA', 'Sample80': 'NA', 'Sample20': 'NA', 'Sample38': 'NA', 'Sample26': 'NA', 'Sample94': 'NA', 'Sample48': 'NA', 'Sample68': 'NA', 'Sample90': 'NA', 'Sample18': 'NA', 'Sample54': 'NA', 'Sample52': 'NA', 'Sample86': 'NA', 'Sample30': 'NA', 'Sample78': 'NA', 'Sample84': 'NA', 'Sample82': 'NA', 'Sample6': 'NA', 'Sample96': 'NA', 'Sample72': 'NA', 'Sample56': 'NA', 'Sample64': 'NA', 'Sample58': 'NA', 'Sample92': 'NA', 'Sample62': 'NA', 'Sample8': 'NA', 'Sample24': 'NA', 'Sample12': 'NA', 'Sample16': 'NA', 'Sample88': 'NA', 'Sample22': 'NA', 'Sample42': 'NA', 'Sample76': 'NA', 'Sample28': 'NA', 'Sample74': 'NA', 'Sample50': 'NA', 'Sample60': 'NA', 'Sample10': 'NA', 'Sample36': 'NA', 'Sample34': 'NA', 'Sample40': 'NA', 'Sample66': 'NA', 'Sample14': 'NA', 'Sample32': 'NA', 'Sample70': 'NA', 'Sample4': 'NA', 'Sample1': '28.00'}
        self.assertEqual(metrics['msi'], expected_msis)

        expected_msi_statuses = {'Sample46': 'NA', 'Sample44': 'NA', 'Sample80': 'NA', 'Sample20': 'NA', 'Sample38': 'NA', 'Sample26': 'NA', 'Sample94': 'NA', 'Sample48': 'NA', 'Sample68': 'NA', 'Sample90': 'NA', 'Sample18': 'NA', 'Sample54': 'NA', 'Sample52': 'NA', 'Sample86': 'NA', 'Sample30': 'NA', 'Sample78': 'NA', 'Sample84': 'NA', 'Sample82': 'NA', 'Sample6': 'NA', 'Sample96': 'NA', 'Sample72': 'NA', 'Sample56': 'NA', 'Sample64': 'NA', 'Sample58': 'NA', 'Sample92': 'NA', 'Sample62': 'NA', 'Sample8': 'NA', 'Sample24': 'NA', 'Sample12': 'NA', 'Sample16': 'NA', 'Sample88': 'NA', 'Sample22': 'NA', 'Sample42': 'NA', 'Sample76': 'NA', 'Sample28': 'NA', 'Sample74': 'NA', 'Sample50': 'NA', 'Sample60': 'NA', 'Sample10': 'NA', 'Sample36': 'NA', 'Sample34': 'NA', 'Sample40': 'NA', 'Sample66': 'NA', 'Sample14': 'NA', 'Sample32': 'NA', 'Sample70': 'NA', 'Sample4': 'NA', 'Sample1': 'Instable'}
        self.assertEqual(metrics['msi_status'], expected_msi_statuses)

    ##################################################################
    #
    #  ######## ##     ## ##       ##           #######
    #  ##       ##     ## ##       ##          ##     ##
    #  ##       ##     ## ##       ##                 ##
    #  ######   ##     ## ##       ##           #######
    #  ##       ##     ## ##       ##          ##
    #  ##       ##     ## ##       ##          ##
    #  ##        #######  ######## ########    #########
    #
    #################################################################
    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_run_worflow_two_mafs(self):
        """
        Test that the workflow works correctly when run with two maf files
        """
        data_clinical_file = os.path.join(self.DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(self.DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        microsatellites_file = self.DATA_SETS['demo']['microsatellites_file']

        tumor_bam1 = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample1.rg.md.abra.printreads.bam")
        normal_bam1 = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample2.rg.md.abra.printreads.bam")

        tumor_bam2 = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample4.rg.md.abra.printreads.bam")
        normal_bam2 = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample3.rg.md.abra.printreads.bam")

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
            "microsatellites_file": {
                "path": microsatellites_file,
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
            ],
            # these must be in the same order as pairs
            "normal_bam_files": [
                {'class': 'File', 'path': normal_bam1},
                {'class': 'File', 'path': normal_bam2}
            ],
            "tumor_bam_files": [
                {'class': 'File', 'path': tumor_bam1},
                {'class': 'File', 'path': tumor_bam2}
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_dir': ODir(name = 'analysis', dir = output_dir, items = [
                OFile(name = 'Proj_08390_G.gene.cna.txt', size = 170207, hash = 'adf367afca81f3bb6e75272c86c38c4ad20246a9'),
                OFile(name = 'Proj_08390_G.muts.maf', size = 63615, hash = '213939883b558f956cbe007909c6d68e60a3add2'),
                OFile(name = 'Proj_08390_G.muts.share.maf', size = 10795, hash = '2dafe49600cac4fbc00e230af45ee2eb258aa6ac'),
                OFile(name = 'Proj_08390_G.seg.cna.txt', size = 2703, hash = '8995133edb5b4547371b0fdea9925848eb480fbc'),
                OFile(name = 'Proj_08390_G.svs.maf', size = 35595, hash = '5c2a63fc01980550108e58079a8b689d53c97d8c')
            ]),
                'facets_dir': ODir(name = 'facets', dir = output_dir, items = [
                    ODir(name = 'Sample1.Sample2', items = [
                        OFile(name = 'Sample1.Sample2_hisens.ccf.portal.maf', size = 19822273, hash = 'fa311e483303fbda7778622af7be57c7d62aac76'),
                        OFile(name = 'Sample1.arm_level.txt', size = 1986, hash = '4f1d6620f3c613e40a147338739d3267e891dda7'),
                        OFile(name = 'Sample1.txt', size = 560, hash = '8abe9bed5d1cb807b451c0647f658394d3f96708'),
                        OFile(name = 'Sample1.gene_level.txt', size = 252793, hash = '6b31bff7758586ec39c94a0c2ddb8201eda5e0ce'),
                        OFile(name = 'Sample1_hisens.cncf.txt', size = 5774, hash = 'a183b43abb9b7bb8f67ec6a0411539643d1275fe'),
                        OFile(name = 'Sample1_hisens.rds', size = 550930, hash = 'd46031ce35efd7341146cac0db809c8cb75c97c2'),
                        OFile(name = 'Sample1_hisens.seg', size = 1786, hash = '1e0f1ad1cc751d53f02ad19717921027c3c96dfb'),
                        OFile(name = 'Sample1_purity.rds', size = 550388, hash = '265c7b94df391ae12e3bc1cc3ec4e763629a4f0a'),
                        OFile(name = 'Sample1_purity.seg', size = 1284, hash = '36904773ac5f194eac7f816a5fa9b055931bb48c'),
                        OFile(name = 'Sample1.qc.txt', size = 1340, hash = '6fdefe96a23c0bc52050c46a16625fca603e4eb6'),
                    ]),
                    ODir(name = 'Sample4.Sample3', items = [
                        OFile(name = 'Sample4.Sample3_hisens.ccf.portal.maf', size = 1136054, hash = '86af92550b3e4252f950223f997579bfe651b65f'),
                        OFile(name = 'Sample4.arm_level.txt', size = 1424, hash = '09b97c42a9875d642dae4d030712264f1b92b550'),
                        OFile(name = 'Sample4.txt', size = 602, hash = '1bb6cf1ea21bc774037913633e5e06772d791b35'),
                        OFile(name = 'Sample4.gene_level.txt', size = 144319, hash = '4a1f580aebf86a3a2ce4d5d27bc1caa3655e8788'),
                        OFile(name = 'Sample4_hisens.cncf.txt', size = 4802, hash = 'a51eac1442aac591412b79305f003d464322140e'),
                        OFile(name = 'Sample4_hisens.rds', size = 186732, hash = '33286048372616b019fcf52587b25bbbf5d1ae08'),
                        OFile(name = 'Sample4_hisens.seg', size = 1584, hash = 'ea32de058625243fdc9f48d9cbab746444549e2e'),
                        OFile(name = 'Sample4_purity.rds', size = 186513, hash = '1eda8dc016cf8322cab1190ea4f723ecb33c4440'),
                        OFile(name = 'Sample4_purity.seg', size = 1276, hash = '6be6775133173d609e176f16a70c5a3c2aeb9027'),
                        OFile(name = 'Sample4.qc.txt', size = 1330, hash = '0c17bd63f43a6726682dec93a58cbcb3ccfd20e1'),
                    ])
                ]),
                'facets_failed_pairs': [],
                'portal_dir': ODir(name = 'portal', dir = output_dir, items = [
                    OFile(name = 'meta_clinical_sample.txt', size = 140, hash = '4c567d81c3b17a76c324fd3e2f73793a6e804f65'),
                    OFile(name = 'data_clinical_patient.txt', size = 643, hash = '9417dcabddd6ab2cbe98167bccd9b9e4fa182562'),
                    OFile(name = 'data_clinical_sample.txt', size = 9161, hash = '519bf38651910dd2954ba959d845962da377f1c0'),
                    OFile(name = 'meta_study.txt', size = 152, hash = '2b0a5fd1a97329adf7c3b1596c84cd6567059a95'),
                    OFile(name = 'meta_clinical_patient.txt', size = 142, hash = '9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f'),
                    OFile(name = 'meta_CNA.txt', size = 270, hash = 'a9bf16f6a0490b19e611e8814b85f7bf1d52417a'),
                    OFile(name = 'meta_fusions.txt', size = 227, hash = '77649e888bafc6a4ed61261d1c46d2f238e1c32b'),
                    OFile(name = 'meta_mutations_extended.txt', size = 253, hash = 'fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d'),
                    OFile(name = 'Proj_08390_G_meta_cna_hg19_seg.txt', size = 200, hash = '59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e'),
                    OFile(name = 'data_CNA.txt', size = 6727, hash = '4795ba3d6c30a028d21bafbb55227ede294b61d5'),
                    OFile(name = 'data_CNA.ascna.txt', size = 8796, hash = '1081a3ab9830e0d1eaddafb45e93956264fdfba9'),
                    OFile(name = 'data_mutations_extended.txt', size = 8382, hash = 'e688086833d015e38e4dea2e207552190bae5d21'),
                    OFile(name = 'Proj_08390_G_data_cna_hg19.seg', size = 2703, hash = '8995133edb5b4547371b0fdea9925848eb480fbc'),
                    OFile(name = 'data_fusions.txt', size = 99, hash = 'c16f763b248813fcdde76f7486f1ddc4e9856038'),
                    ODir(name = 'case_lists', items = [
                        OFile(name = 'cases_all.txt', size = 616, hash = 'b9e43289cec5603b0886b5e8507c8d019387c125'),
                        OFile(name = 'cases_cnaseq.txt', size = 696, hash = 'b87e2da8dce0fddbadec348efe2986519b2a794b'),
                        OFile(name = 'cases_cna.txt', size = 628, hash = '053481a8299e9430117f8e45e081aa7ec21033a6'),
                        OFile(name = 'cases_sequenced.txt', size = 641, hash = 'ef9f5aef03c2527bf576470168660557ca1c7cc9'),
                    ]),
                    OFile(name = 'report.html')
                ])
            }

        for item in output_json['portal_dir']['listing']:
            if item['basename'] == 'report.html':
                item.pop('checksum')
                item.pop('size')

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
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tgenome_doubled\tASCN_WGD\tCMO_TMB_SCORE\tCMO_MSI_SCORE\tCMO_MSI_STATUS\n',
        '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tSAMPLE_COVERAGE\tPROJECT_PI\tREQUEST_PI\tASCN_PURITY\tASCN_PLOIDY\tASCN_VERSION\tgenome_doubled\tASCN_WGD\tCMO_TMB_SCORE\tCMO_MSI_SCORE\tCMO_MSI_STATUS\n',
        '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\tSTRING\tSTRING\tNUMBER\tNUMBER\tSTRING\n',
        '#1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t0\t0\n'
        ]
        self.assertEqual(comments, expected_comments)

        metrics = {}
        metrics['tmb'] = {}
        metrics['msi'] = {}
        metrics['msi_status'] = {}
        for record in records:
            sample_id = record['SAMPLE_ID']
            metrics['tmb'][sample_id] = record['CMO_TMB_SCORE']
            metrics['msi'][sample_id] = record['MSI_SCORE']
            metrics['msi_status'][sample_id] = record['MSI_STATUS']

        expected_tmbs = {
        'Sample46': 'NA', 'Sample44': 'NA', 'Sample80': 'NA', 'Sample20': 'NA', 'Sample38': 'NA', 'Sample26': 'NA', 'Sample94': 'NA', 'Sample48': 'NA', 'Sample68': 'NA', 'Sample90': 'NA', 'Sample18': 'NA', 'Sample54': 'NA', 'Sample52': 'NA', 'Sample86': 'NA', 'Sample30': 'NA', 'Sample78': 'NA', 'Sample84': 'NA', 'Sample82': 'NA', 'Sample6': 'NA', 'Sample96': 'NA', 'Sample72': 'NA', 'Sample56': 'NA', 'Sample64': 'NA', 'Sample58': 'NA', 'Sample92': 'NA', 'Sample62': 'NA', 'Sample8': 'NA', 'Sample24': 'NA', 'Sample12': 'NA', 'Sample16': 'NA', 'Sample88': 'NA', 'Sample22': 'NA', 'Sample42': 'NA', 'Sample76': 'NA', 'Sample28': 'NA', 'Sample74': 'NA', 'Sample50': 'NA', 'Sample60': 'NA', 'Sample10': 'NA', 'Sample36': 'NA', 'Sample34': 'NA', 'Sample40': 'NA', 'Sample66': 'NA', 'Sample14': 'NA', 'Sample32': 'NA', 'Sample70': 'NA', 'Sample4': '5.5', 'Sample1': '47.5'
        }
        self.assertEqual(metrics['tmb'], expected_tmbs)

        expected_msis = {'Sample46': 'NA', 'Sample44': 'NA', 'Sample80': 'NA', 'Sample20': 'NA', 'Sample38': 'NA', 'Sample26': 'NA', 'Sample94': 'NA', 'Sample48': 'NA', 'Sample68': 'NA', 'Sample90': 'NA', 'Sample18': 'NA', 'Sample54': 'NA', 'Sample52': 'NA', 'Sample86': 'NA', 'Sample30': 'NA', 'Sample78': 'NA', 'Sample84': 'NA', 'Sample82': 'NA', 'Sample6': 'NA', 'Sample96': 'NA', 'Sample72': 'NA', 'Sample56': 'NA', 'Sample64': 'NA', 'Sample58': 'NA', 'Sample92': 'NA', 'Sample62': 'NA', 'Sample8': 'NA', 'Sample24': 'NA', 'Sample12': 'NA', 'Sample16': 'NA', 'Sample88': 'NA', 'Sample22': 'NA', 'Sample42': 'NA', 'Sample76': 'NA', 'Sample28': 'NA', 'Sample74': 'NA', 'Sample50': 'NA', 'Sample60': 'NA', 'Sample10': 'NA', 'Sample36': 'NA', 'Sample34': 'NA', 'Sample40': 'NA', 'Sample66': 'NA', 'Sample14': 'NA', 'Sample32': 'NA', 'Sample70': 'NA', 'Sample4': '11.76', 'Sample1': '28.00'}
        self.assertEqual(metrics['msi'], expected_msis)

        expected_msi_statuses = {'Sample46': 'NA', 'Sample44': 'NA', 'Sample80': 'NA', 'Sample20': 'NA', 'Sample38': 'NA', 'Sample26': 'NA', 'Sample94': 'NA', 'Sample48': 'NA', 'Sample68': 'NA', 'Sample90': 'NA', 'Sample18': 'NA', 'Sample54': 'NA', 'Sample52': 'NA', 'Sample86': 'NA', 'Sample30': 'NA', 'Sample78': 'NA', 'Sample84': 'NA', 'Sample82': 'NA', 'Sample6': 'NA', 'Sample96': 'NA', 'Sample72': 'NA', 'Sample56': 'NA', 'Sample64': 'NA', 'Sample58': 'NA', 'Sample92': 'NA', 'Sample62': 'NA', 'Sample8': 'NA', 'Sample24': 'NA', 'Sample12': 'NA', 'Sample16': 'NA', 'Sample88': 'NA', 'Sample22': 'NA', 'Sample42': 'NA', 'Sample76': 'NA', 'Sample28': 'NA', 'Sample74': 'NA', 'Sample50': 'NA', 'Sample60': 'NA', 'Sample10': 'NA', 'Sample36': 'NA', 'Sample34': 'NA', 'Sample40': 'NA', 'Sample66': 'NA', 'Sample14': 'NA', 'Sample32': 'NA', 'Sample70': 'NA', 'Sample4': 'Instable', 'Sample1': 'Instable'}
        self.assertEqual(metrics['msi_status'], expected_msi_statuses)

if __name__ == "__main__":
    unittest.main()
