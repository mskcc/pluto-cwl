#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import os
import unittest

if __name__ != "__main__":
    from .tools import TmpDirTestCase, run_cwl
    from .settings import CWL_DIR

if __name__ == "__main__":
    from tools import TmpDirTestCase, run_cwl
    from settings import CWL_DIR

cwl_file = os.path.join(CWL_DIR, 'add_header.cwl')

class TestAddHeader(TmpDirTestCase):
    def test_add_header(self):
        """
        Test case for adding a header to a file
        """
        self.maxDiff = None
        
        input_file = os.path.join(self.tmpdir, "input.txt")
        with open(input_file, "w") as f:
            f.write("foo")
        header_str = "HEADER"

        input_json = {
            "input_file": {
                  "class": "File",
                  "path": input_file
                },
            "header_str":  header_str,
            }

        output_json, output_dir = run_cwl(
            testcase = self,
            tmpdir = self.tmpdir,
            input_json = input_json,
            cwl_file = cwl_file,
            print_command = False,
            )

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.txt'),
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$01838a0977d542fb12680e271393e1d4baaefa8f',
                'size': 10,
                'path':  os.path.join(output_dir,'output.txt')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as f:
            lines = [ l.strip() for l in f ]
        expected_lines = ['HEADER', 'foo']
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
