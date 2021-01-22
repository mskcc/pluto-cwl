#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run all the unit tests

NOTE: This does not seem to work correctly, it skips some tests located in the `tests/` directory, 
so do not use this and instead use the `make test` recipes from the Makefile
TODO: figure out what is wrong with this script
"""

import unittest
import sys

if __name__ == "__main__":
    loader = unittest.TestLoader()
    start_dir = '.'
    suite = loader.discover(start_dir)

    runner = unittest.TextTestRunner()
    ret = not runner.run(suite).wasSuccessful()
    sys.exit(ret)
