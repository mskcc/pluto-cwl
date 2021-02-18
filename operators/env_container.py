#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the example workflow
"""
from pluto.tools import CWLFile
from .classes import Operator

class EnvContainerCWL(Operator):
    cwl_file = CWLFile('env_container.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.input = {}
