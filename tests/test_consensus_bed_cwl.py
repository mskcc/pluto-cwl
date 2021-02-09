#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the consensus bed cwl
"""
import os
import sys
import unittest
from collections import OrderedDict

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
sys.path.pop(0)

class TestConsensusBed(PlutoTestCase):
    cwl_file = CWLFile('consensus_bed.cwl')

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
        ('End_Position', '62321135')
        ])
        self.maf_row2 = OrderedDict([
        ('Hugo_Symbol', 'FAM46C'),
        ('Entrez_Gene_Id', '54855'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '1'),
        ('Start_Position', '118166398'),
        ('End_Position', '118166398')
        ])
        self.maf_row3 = OrderedDict([
        ('Hugo_Symbol', 'IL7R'),
        ('Entrez_Gene_Id', '3575'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '5'),
        ('Start_Position', '35876484'),
        ('End_Position', '35876484')
        ])
        self.maf_row4 = OrderedDict([
        ('Hugo_Symbol', 'KMT2C'),
        ('Entrez_Gene_Id', '58508'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '7'),
        ('Start_Position', '151845367'),
        ('End_Position', '151845367')
        ])
        self.maf_row5 = OrderedDict([
        ('Hugo_Symbol', 'MET'),
        ('Entrez_Gene_Id', '4233'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '7'),
        ('Start_Position', '116418998'),
        ('End_Position', '116418998')
        ])
        self.maf_row6 = OrderedDict([
        ('Hugo_Symbol', 'MAP2K4'),
        ('Entrez_Gene_Id', '6416'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '17'),
        ('Start_Position', '11998935'),
        ('End_Position', '11998935')
        ])

        maf_rows1 = [ self.maf_row1, self.maf_row2, self.maf_row3, self.maf_row4 ]
        maf_rows2 = [ self.maf_row1, self.maf_row2, self.maf_row5, self.maf_row6 ]
        maf_lines1 = []
        maf_lines2 = []
        for comment in self.comments:
            maf_lines1.append(comment[0] + '\n')
            maf_lines2.append(comment[0] + '\n')
        header = '\t'.join([ k for k in maf_rows1[0].keys() ])
        header += '\n'
        maf_lines1.append(header)
        maf_lines2.append(header)
        for row in maf_rows1:
            values = [ v for v in row.values() ]
            line = '\t'.join(values)
            line += '\n'
            maf_lines1.append(line)
        for row in maf_rows2:
            values = [ v for v in row.values() ]
            line = '\t'.join(values)
            line += '\n'
            maf_lines2.append(line)
        self.maf1 = os.path.join(self.tmpdir, "input1.maf")
        self.maf2 = os.path.join(self.tmpdir, "input2.maf")
        with open(self.maf1, "w") as fout:
            for line in maf_lines1:
                fout.write(line)
        with open(self.maf2, "w") as fout:
            for line in maf_lines2:
                fout.write(line)

    def test_consensus_bed_workflow(self):
        """
        """
        self.maxDiff = None
        self.input = {
            'maf_files': [
                {'class': 'File', 'path': self.maf1},
                {'class': 'File', 'path': self.maf2},
            ]
        }
        output_json, output_dir = self.run_cwl()
        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'merged.bed'),
                'basename': 'merged.bed',
                'class': 'File',
                'checksum': 'sha1$f831dc91b70c02f10f69da2dae21de57d580b654',
                'size': 149,
                'path':  os.path.join(output_dir,'merged.bed')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = output_json['output_file']['path']

        with open(output_file) as f:
            lines = [ l for l in f ]

        expected_lines = [
            '1\t118166398\t118166398\n',
            '5\t35876484\t35876484\n',
            '7\t116418998\t116418998\n',
            '17\t11998935\t11998935\n',
            '20\t62321135\t62321135\n',
            '7\t151845367\t151845367\n',
            '20\t62321135\t62321135\n'
        ]
        self.assertDictEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
