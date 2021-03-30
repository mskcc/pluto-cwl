#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the msi CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class MSI(Operator):
    cwl_file = CWLFile('msi.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input = generate_input(
            self.args,
            File_keys = ['n', 't', 'd'],
            File_keys_abspath = True
        )
