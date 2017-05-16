#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

from __future__ import print_function
from sys import exit

from exception import TestFailed

class TestEngine:
    @staticmethod
    def run( testcase ):
        try:
            success = testcase.performTest()
            print( '[PASS] {message}'.format( message=success ) )
        except TestFailed as ex:
            exit( '[FAIL] {message}'.format( message=ex ) )
