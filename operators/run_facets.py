#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the run-facets-wrapper CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class RunFacets(Operator):
    cwl_file = CWLFile('run-facets-wrapper.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input = generate_input(
            self.args,
            File_keys = ['snp_pileup'],
            File_keys_abspath = True
        )
