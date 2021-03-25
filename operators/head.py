#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the head cwl
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class HeadCWL(Operator):
    cwl_file = CWLFile('head.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.input = {}

        self.input = generate_input(
            self.args,
            File_keys = ['input_file'], File_keys_abspath = True
        )
