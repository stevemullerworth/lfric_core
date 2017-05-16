#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

from abc import ABCMeta, abstractmethod
import subprocess
import sys

class Test():
  __metaclass__ = ABCMeta

  def __init__( self, command=sys.argv[1] ):
    if type(command) is not list:
      command = [command]

    self.process = subprocess.Popen( command,
                                     stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE )

  def performTest( self ):
    return self.test( self.process )

  @abstractmethod
  def test( self, process ):
    pass

class MpiTest(Test):
  __metaclass__ = ABCMeta

  def __init__( self, command=sys.argv[1], cores=4 ):
    mpiCommand = ['mpiexec', '-n', str(cores)]
    if type(command) is list:
      mpiCommand.extend( command )
    else:
      mpiCommand.append( command )

    super(MpiTest, self).__init__( mpiCommand )
