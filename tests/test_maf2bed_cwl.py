#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the maf2bed cwl
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

class TestMaf2Bed(PlutoTestCase):
    cwl_file = CWLFile('maf2bed.cwl')

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

    def test_maf2bed(self):
        maf_rows = [ self.maf_row1, self.maf_row2, self.maf_row3, self.maf_row4, self.maf_row5, self.maf_row6 ]
        maf_lines = []
        for comment in self.comments:
            maf_lines.append(comment[0] + '\n')
        header = '\t'.join([ k for k in maf_rows[0].keys() ])
        header += '\n'
        maf_lines.append(header)
        for row in maf_rows:
            values = [ v for v in row.values() ]
            line = '\t'.join(values)
            line += '\n'
            maf_lines.append(line)
        maf = os.path.join(self.tmpdir, "input1.maf")
        with open(maf, "w") as fout:
            for line in maf_lines:
                fout.write(line)

        self.input = {
            "maf_file": { "class": "File", "path": maf }
            # "output_filename": "output.bed"
        }
        output_json, output_dir = self.run_cwl()
        path = output_json['output_file'].pop('path')
        location = output_json['output_file'].pop('location')
        basename = output_json['output_file'].pop('basename')

        self.assertTrue(basename.startswith('_maf2bed_merged'))

        expected_output = {
            'output_file': {
                # 'location': 'file://' + os.path.join(output_dir,'output.bed'),
                # 'basename': 'output.bed',
                'class': 'File',
                'checksum': 'sha1$3fdf885d0f863d366c3a03154e358c64ea3465a3',
                'size': 128,
                # 'path': os.path.join(output_dir,'output.bed'),
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        with open(path) as fin:
            lines = [ line for line in fin ]

        expected_lines = [
            '1\t118166398\t118166398\n',
            '5\t35876484\t35876484\n',
            '7\t116418998\t116418998\n',
            '7\t151845367\t151845367\n',
            '17\t11998935\t11998935\n',
            '20\t62321135\t62321135\n'
        ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
