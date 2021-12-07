#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Integration tests for the workflow_with_facets.cwl using medium sized dataset

Usage:

    $ INTEGRATION_TESTS=True USE_LSF=True CWL_ENGINE=toil python tests/test_workflow_with_facets.medium.py

"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase
from pluto.settings import ENABLE_INTEGRATION_TESTS
# from pluto.serializer import OFile, ODir
sys.path.pop(0)

from fixtures import WORKFLOW_MEDIUM_JSON

class TestWorkflowWithFacetsMedium(PlutoTestCase):
    cwl_file = 'workflow_with_facets.cwl'

    @unittest.skipIf(ENABLE_INTEGRATION_TESTS!=True, "is a large integration test")
    def test_medium_workflow(self):
        """
        $ PRESERVE_TEST_DIR=True USE_LSF=True PRINT_COMMAND=True CWL_ENGINE=toil python tests/test_workflow_with_facets.py TestWorkflowWithFacets.test_medium_workflow
        ---------
        Ran 1 test in 8698.861s
        """
        input_json = WORKFLOW_MEDIUM_JSON
        output_json, output_dir = self.run_cwl(input = input_json, input_is_file = True)

if __name__ == "__main__":
    unittest.main()
