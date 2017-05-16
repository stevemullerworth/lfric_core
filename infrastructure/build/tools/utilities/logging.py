#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
# A library of logging tools.

from __future__ import print_function

import abc
import sys

###############################################################################
# Abstract logging provider class
class Logger():
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def logEvent( self, description ):
        raise NotImplementedError( 'Logger.logEvent not implemented' )

###############################################################################
# Silent (null) logger
class NoLogger( Logger ):
    def logEvent( self, description ):
        pass

###############################################################################
# Print to file logger
class PrintLogger( Logger ):
    def __init__(self, stream ):
        self.stream = stream

    def logEvent( self, description ):
        print( '{}'.format( description ), file=self.stream )
        sys.stdout.flush()
