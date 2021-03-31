#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the samples_fillout_index_workflow cwl
"""
import os
import sys
import unittest
from collections import OrderedDict
import shutil

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile, TableReader
sys.path.pop(0)

class TestSamplesFilloutIndex(PlutoTestCase):
    cwl_file = CWLFile('samples_fillout_index_workflow.cwl')

    def setUp(self):
        super().setUp()
        self.comments = [
        ['# comment 1'],
        ['# comment 2']
        ]
        self.maf_row1 = OrderedDict([
        ('Hugo_Symbol', 'RTEL1'),
        ('Entrez_Gene_Id', '51750'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '20'),
        ('Start_Position', '62321135'),
        ('End_Position', '62321135'),
        ('Variant_Classification', 'Silent'), # all the columns after this do not matter
        ('Reference_Allele', 'G'),
        ('Tumor_Seq_Allele1', 'G'),
        ('Tumor_Seq_Allele2', 'A'),
        ('n_alt_count', '1'),
        ('Matched_Norm_Sample_Barcode', '.'),
        ('t_alt_count', '142'),
        ('t_ref_count', '511'),
        ('n_ref_count', '212'),
        ('Tumor_Sample_Barcode', '.')
        ])
        self.maf_row2 = OrderedDict([
        ('Hugo_Symbol', 'FAM46C'),
        ('Entrez_Gene_Id', '54855'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '1'),
        ('Start_Position', '118166398'),
        ('End_Position', '118166398'),
        ('Variant_Classification', 'Silent'),
        ('Reference_Allele', 'G'),
        ('Tumor_Seq_Allele1', 'G'),
        ('Tumor_Seq_Allele2', 'A'),
        ('n_alt_count', '1'),
        ('Matched_Norm_Sample_Barcode', '.'),
        ('t_alt_count', '142'),
        ('t_ref_count', '511'),
        ('n_ref_count', '212'),
        ('Tumor_Sample_Barcode', '.')
        ])
        self.maf_row3 = OrderedDict([
        ('Hugo_Symbol', 'IL7R'),
        ('Entrez_Gene_Id', '3575'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '5'),
        ('Start_Position', '35876484'),
        ('End_Position', '35876484'),
        ('Variant_Classification', 'Silent'),
        ('Reference_Allele', 'G'),
        ('Tumor_Seq_Allele1', 'G'),
        ('Tumor_Seq_Allele2', 'A'),
        ('n_alt_count', '1'),
        ('Matched_Norm_Sample_Barcode', '.'),
        ('t_alt_count', '142'),
        ('t_ref_count', '511'),
        ('n_ref_count', '212'),
        ('Tumor_Sample_Barcode', '.')
        ])
        self.maf_row4 = OrderedDict([
        ('Hugo_Symbol', 'KMT2C'),
        ('Entrez_Gene_Id', '58508'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '7'),
        ('Start_Position', '151845367'),
        ('End_Position', '151845367'),
        ('Variant_Classification', 'Silent'),
        ('Reference_Allele', 'G'),
        ('Tumor_Seq_Allele1', 'G'),
        ('Tumor_Seq_Allele2', 'A'),
        ('n_alt_count', '1'),
        ('Matched_Norm_Sample_Barcode', '.'),
        ('t_alt_count', '142'),
        ('t_ref_count', '511'),
        ('n_ref_count', '212'),
        ('Tumor_Sample_Barcode', '.')
        ])

        rows1 = [ self.maf_row1, self.maf_row2 ]
        lines1 = self.dicts2lines(rows1, comment_list = self.comments)
        self.maf1 = self.write_table(tmpdir = self.tmpdir, filename = "1.maf", lines = lines1)

        rows2 = [ self.maf_row3, self.maf_row4 ]
        lines2 = self.dicts2lines(rows2, comment_list = self.comments)
        self.maf2 = self.write_table(tmpdir = self.tmpdir, filename = "2.maf", lines = lines2)

        rows3 = [ self.maf_row1, self.maf_row3 ]
        lines3 = self.dicts2lines(rows3, comment_list = self.comments)
        self.maf3 = self.write_table(tmpdir = self.tmpdir, filename = "3.maf", lines = lines3)

        # rows4 = [ self.maf_row2, self.maf_row4 ]
        # lines4 = self.dicts2lines(rows4, comment_list = self.comments)
        # self.maf4 = self.write_table(tmpdir = self.tmpdir, filename = "4.maf", lines = lines4)

        # copy some bam files so they will not have .bai files present
        self.bam1 = os.path.join(self.tmpdir, "1.bam")
        shutil.copyfile(os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam"), self.bam1)

        self.bam2 = os.path.join(self.tmpdir, "2.bam")
        shutil.copyfile(os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam"), self.bam2)

        # this one still has its .bai file # Sample36.Sample35
        self.bam3 = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample35.rg.md.abra.printreads.bam")


    def test_run_fillout_workflow(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        """
        self.maxDiff = None

        self.input = {
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
            "sample_ids": ["Sample35"],
            "bam_files": [
                { "class": "File", "path": self.bam3 }
            ],
            "maf_files": [
                { "class": "File", "path": self.maf3 }
            ],
            "unindexed_sample_ids": ["Sample24", "Sample23"],
            "unindexed_bam_files": [
                { "class": "File", "path": self.bam1 },
                { "class": "File", "path": self.bam2 }
            ],
            "unindexed_maf_files": [
                { "class": "File", "path": self.maf1 },
                { "class": "File", "path": self.maf2 }
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.maf'),
                'basename': 'output.maf',
                'class': 'File',
                'checksum': 'sha1$9c692bdffe67c88af06c528ad82ee25ce95a633c',
                'size': 2447,
                'path':  os.path.join(output_dir,'output.maf')
                }
            }
        self.assertEqual(output_json, expected_output)

        output_file = output_json['output_file']['path']
        reader = TableReader(output_file)
        comments = reader.comment_lines
        fieldnames = reader.get_fieldnames()
        records = [ rec for rec in reader.read() ]

        expected_records =  [
        {'Hugo_Symbol': 'FAM46C', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '1', 'Start_Position': '118166398', 'End_Position': '118166398', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample35', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '0', 't_alt_count': '9', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '9', 't_variant_frequency': '1', 't_total_count_forward': '4', 't_ref_count_forward': '0', 't_alt_count_forward': '4'},
        {'Hugo_Symbol': 'IL7R', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '5', 'Start_Position': '35876484', 'End_Position': '35876484', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample35', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '0', 't_alt_count': '0', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '21', 't_variant_frequency': '0', 't_total_count_forward': '9', 't_ref_count_forward': '0', 't_alt_count_forward': '0'},
        {'Hugo_Symbol': 'KMT2C', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '7', 'Start_Position': '151845367', 'End_Position': '151845367', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample35', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '10', 't_alt_count': '0', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '10', 't_variant_frequency': '0', 't_total_count_forward': '5', 't_ref_count_forward': '5', 't_alt_count_forward': '0'},
        {'Hugo_Symbol': 'RTEL1', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '20', 'Start_Position': '62321135', 'End_Position': '62321135', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample35', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '20', 't_alt_count': '0', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '20', 't_variant_frequency': '0', 't_total_count_forward': '10', 't_ref_count_forward': '10', 't_alt_count_forward': '0'},
        {'Hugo_Symbol': 'FAM46C', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '1', 'Start_Position': '118166398', 'End_Position': '118166398', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample24', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '0', 't_alt_count': '41', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '41', 't_variant_frequency': '1', 't_total_count_forward': '22', 't_ref_count_forward': '0', 't_alt_count_forward': '22'},
        {'Hugo_Symbol': 'IL7R', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '5', 'Start_Position': '35876484', 'End_Position': '35876484', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample24', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '0', 't_alt_count': '0', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '52', 't_variant_frequency': '0', 't_total_count_forward': '26', 't_ref_count_forward': '0', 't_alt_count_forward': '0'},
        {'Hugo_Symbol': 'KMT2C', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '7', 'Start_Position': '151845367', 'End_Position': '151845367', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample24', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '68', 't_alt_count': '4', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '72', 't_variant_frequency': '0.0555556', 't_total_count_forward': '34', 't_ref_count_forward': '32', 't_alt_count_forward': '2'},
        {'Hugo_Symbol': 'RTEL1', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '20', 'Start_Position': '62321135', 'End_Position': '62321135', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample24', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '129', 't_alt_count': '0', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '129', 't_variant_frequency': '0', 't_total_count_forward': '66', 't_ref_count_forward': '66', 't_alt_count_forward': '0'},
        {'Hugo_Symbol': 'FAM46C', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '1', 'Start_Position': '118166398', 'End_Position': '118166398', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample23', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '0', 't_alt_count': '40', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '40', 't_variant_frequency': '1', 't_total_count_forward': '20', 't_ref_count_forward': '0', 't_alt_count_forward': '20'},
        {'Hugo_Symbol': 'IL7R', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '5', 'Start_Position': '35876484', 'End_Position': '35876484', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample23', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '0', 't_alt_count': '0', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '120', 't_variant_frequency': '0', 't_total_count_forward': '58', 't_ref_count_forward': '0', 't_alt_count_forward': '0'},
        {'Hugo_Symbol': 'KMT2C', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '7', 'Start_Position': '151845367', 'End_Position': '151845367', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample23', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '91', 't_alt_count': '0', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '91', 't_variant_frequency': '0', 't_total_count_forward': '43', 't_ref_count_forward': '43', 't_alt_count_forward': '0'},
        {'Hugo_Symbol': 'RTEL1', 'Entrez_Gene_Id': '.', 'Center': 'mskcc.org', 'NCBI_Build': 'GRCh37', 'Chromosome': '20', 'Start_Position': '62321135', 'End_Position': '62321135', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'A', 'Tumor_Seq_Allele2': '.', 'dbSNP_RS': '.', 'dbSNP_Val_Status': '.', 'Tumor_Sample_Barcode': 'Sample23', 'Matched_Norm_Sample_Barcode': 'Normal', 'Match_Norm_Seq_Allele1': '.', 'Match_Norm_Seq_Allele2': '.', 'Tumor_Validation_Allele1': '.', 'Tumor_Validation_Allele2': '.', 'Match_Norm_Validation_Allele1': '.', 'Match_Norm_Validation_Allele2': '.', 'Verification_Status': '.', 'Validation_Status': '.', 'Mutation_Status': 'UNPAIRED', 'Sequencing_Phase': '.', 'Sequence_Source': '.', 'Validation_Method': '.', 'Score': '.', 'BAM_File': '.', 'Sequencer': '.', 't_ref_count': '184', 't_alt_count': '0', 'n_ref_count': '.', 'n_alt_count': '.', 'Caller': '.', 't_total_count': '184', 't_variant_frequency': '0', 't_total_count_forward': '89', 't_ref_count_forward': '89', 't_alt_count_forward': '0'}
        ]

        self.assertEqual(records, expected_records)

if __name__ == "__main__":
    unittest.main()
