#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
# A library of filename tools.

from __future__ import print_function

import os.path

###############################################################################
def replaceExtension( filename, extension ):
    (base, rubbish) = os.path.splitext( filename )
    return '{}.{}'.format( base, extension)