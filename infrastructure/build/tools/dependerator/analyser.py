#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
# Generate a make file snippet holding build information about a source file.
#
# This snippet consists of a dependency list built up from which modules a
# file "use"s.
#
# Additionally a separate snippet file is created for each "program" unit found
# listing all the modules which go to make that program.
#
# These snippets may then be "include"ed into other make files.

from __future__ import print_function;

from abc import ABCMeta, abstractmethod
import os
import os.path
import re
import shutil
import subprocess

'''
Examine Fortran source and build dependency information for use by "make".
'''

###############################################################################
# Interface for analysers
#
class Analyser(object):
    __metaclass__ = ABCMeta

    ###########################################################################
    # Examine a source file and store dependency information in the database.
    #
    # Arguments:
    #   sourceFilename - The name of the object to scan.
    #
    @abstractmethod
    def analyse( self, sourceFilename ):
        pass

###############################################################################
# Examine a namelist description file for dependencies.
#
class NamelistDescriptionAnalyser(Analyser):
    ###########################################################################
    # Constructor
    # Arguments:
    #   logger   - A logger object to write output to.
    #   database - FileDependencies object to hold details.
    #
    def __init__( self, logger, database ):
        self._logger   = logger
        self._database = database

    ###########################################################################
    # Scan a namelist description file and harvest dependency information.
    #
    # Arguments:
    #   sourceFilename - File object to scan.
    #
    def analyse( self, sourceFilename ):
        if not sourceFilename.endswith( '.nld' ):
            raise Exception( 'File doesn''t look like a namelist description' \
                             + ' file: ' + sourceFilename )

        self._logger.logEvent( '  Scanning ' + sourceFilename )
        self._database.removeFile( sourceFilename )

        with open( sourceFilename, 'r' ) as sourceFile:
            for line in sourceFile:
                match = re.match( r'^\s*namelist\s+(\S+)', line, \
                                  flags=re.IGNORECASE )
                if match is not None:
                    fortranFilename = '{}_configuration_mod.f90' \
                                      .format( match.group( 1 ) )
                    self._database.addFileDependency( fortranFilename, \
                                                      sourceFilename )

###############################################################################
# Examine Fortran source code for module dependencies.
#
class FortranAnalyser(Analyser):
  ###########################################################################
  # Constructor
  #
  # Arguments:
  #   logger        - A Logger object to write output to.
  #   ignoreModules - A list of module names to ignore.
  #   database      - FortranDatabase object to hold details.
  #
  def __init__( self, logger, ignoreModules, database ):
    self._logger          = logger
    self._ignoreModules   = map(str.lower, ignoreModules)
    self._database        = database

    self._ignoreModules.append( 'iso_fortran_env' )
    self._ignoreModules.append( 'omp_lib' )

    self._fpp = os.getenv( 'FPP', None )
    if not self._fpp:
        raise Exception( 'No Fortran preprocessor provided in $FPP' )
    self._fpp = self._fpp.split()

    # Patterns to recognise scoping units
    #
    self._programPattern = re.compile( r'^\s*PROGRAM\s+(\S+)',
                                       flags=re.IGNORECASE )
    self._modulePattern \
  = re.compile( r'^\s*MODULE\s+(?!(?:PROCEDURE|SUBROUTINE|FUNCTION)\s+)(\S+)',
                flags=re.IGNORECASE )
    self._submodulePattern \
            = re.compile( r'^\s*SUBMODULE\s+\((?:([^:]+):)?([^)]+)\)\s+(\S+)',
                          flags=re.IGNORECASE )
    self._subroutinePattern \
                     = re.compile( r'^\s*(MODULE\s+)?SUBROUTINE\s+([^\(\s]+)',
                                   flags=re.IGNORECASE )
    self._functionPattern \
= re.compile( r'^\s*(?:TYPE\s*\(\s*\S+?\s*\)\s*|\S+\s*)?(MODULE\s+)?FUNCTION\s+([^\(\s]+)',
              flags=re.IGNORECASE )
    self._endPattern = re.compile( r'^\s*END(?:\s+(\S+)(?:\s+(\S+))?)?',
                                   flags=re.IGNORECASE )

    # Patterns to recognise dependencies
    #
    self._usePattern = re.compile( r'^\s*USE\s+([^,\s]+)',
                                   flags=re.IGNORECASE )
    self._externalPattern \
                 = re.compile( r'^\s*EXTERNAL\s+([^,\s]+(:?\s*,\s*[^,\s]+)*)',
                               flags=re.IGNORECASE )
    self._pFUnitPattern = re.compile( r'^\s*#\s+\d+\s+".*testSuites.inc"',
                                      flags=re.IGNORECASE )
    self._suitePattern  = re.compile( r'^\s*ADD_TEST_SUITE\(\s*(\S+)\)' )
    self._dependsPattern = re.compile( r'!\s*DEPENDS ON:\s*([^.\s]+)(.o)?',
                                       flags=re.IGNORECASE)

  ###########################################################################
  # Scan a Fortran source file and harvest dependency information.
  #
  # Arguments:
  #   sourceFilename(str) - The name of the object to scan.
  #   preprocessMacros(dict) - Macro name is the key. Value may be None for
  #                            empty macros.
  #
  def analyse( self, sourceFilename, preprocessMacros=None ):
    # Perform any necessary preprocessing
    #
    if sourceFilename.endswith( '.F90' ):
      self._logger.logEvent( '  Preprocessing ' + sourceFilename )
      preprocessCommand = self._fpp
      if preprocessMacros:
        for name, macro in preprocessMacros.iteritems():
          if macro:
            preprocessCommand.append( '-D{}={}'.format( name, macro ) )
          else:
            preprocessCommand.append( '-D{}'.format( name ) )
      preprocessCommand.append( sourceFilename )
      preprocessor = subprocess.Popen( preprocessCommand, \
                                       stdout=subprocess.PIPE, \
                                       stderr=subprocess.PIPE )
      processed_source, errors = preprocessor.communicate()
      if preprocessor.returncode:
        raise subprocess.CalledProcessError( preprocessor.returncode,
                                             " ".join( preprocessCommand ) )
    elif sourceFilename.endswith( '.f90' ):
      with open( sourceFilename, 'r' ) as sourceFile:
        processed_source = sourceFile.read()
    else:
      message = 'File doesn''t look like a Fortran file: {}'
      raise Exception( message.format( sourceFilename ) )

    # A local function to hide the mesiness of adding a dependency to the
    # database.
    #
    def add_dependency( program_unit, prerequisite_unit, reverseLink=False ):
      prerequisite_unit = prerequisite_unit.lower()
      if prerequisite_unit in self._ignoreModules:
        self._logger.logEvent( '      - Ignored 3rd party prerequisite' )
      elif prerequisite_unit in modules:
        self._logger.logEvent( '      - Ignored self' )
      elif prerequisite_unit in dependencies:
        self._logger.logEvent( '      - Ignore duplicate prerequisite' )
      else:
        dependencies.append( prerequisite_unit )
        self._database.addModuleCompileDependency( program_unit, \
                                                   prerequisite_unit )
        if reverseLink:
          self._database.addModuleLinkDependency( prerequisite_unit, \
                                                  program_unit )
        else:
          self._database.addModuleLinkDependency( program_unit, \
                                                  prerequisite_unit )

    # Read lines from the file, concatenate at continuation markers and
    # split comments off.
    #
    # TODO: This is complex. That complexity comes from the need
    #       to preserve comments. This is needed to support "depends on"
    #       comments. Ergo, once "depends on" is gone we can ignore comments
    #       and this becomes a lot simpler.
    #
    def lines_of_code( source ):
      line_number = 0
      code    = ''
      comment = ''
      continuing = False
      for line in source.splitlines(): # Loop over every line in the source
        line_number += 1
        state = 'indent' # Each line starts in the "indent" state
        index = -1
        code_start = 0
        code_end = len(line)
        comment_start = len(line)
        continuation = False
        for character in line: # Scan every character in a line
          index += 1
          if state == 'indent': ##############################################
            if character == '&': # Start of continuation line
              if not continuing:
                message = 'Found continuation marker at start of line when there was none ending previous line'
                raise Exception( message )
            elif character == '!': # Line contains only a comment
              comment_start = index
              state = 'comment'
            elif character != ' ': # Start of code located
              code_start = index
              state = 'code'
          elif state == 'code': ##############################################
            if character == '"': # String opened with double quote
              state == 'double'
            elif character == "'": # String opened with single quote
              state == 'single'
            elif character == '&': # Line continues on next line
              code_end = index
              continuing = True
              continuation = True
              state = 'continue'
            elif character == '!': # The remainder of the line is a comment
              comment_start = index
              state = 'comment'
          elif state == 'double': ############################################
            if character == '"': # Quoted string has ended
              state = 'code'
          elif state == 'single': ############################################
            if character == "'": # Quoted string has ended
              state = 'code'
          elif state == 'continue': ##########################################
            if character == '!': # There is a comment after the continuation
              state = 'comment'
          elif state == 'comment': ###########################################
            if index-1 < code_end: # If we have not already ended the code
              code_end = index-1   # Mark it as ended.
            break

        code    += ' ' + line[code_start:code_end]
        comment += ' ' + line[comment_start:]

        if line[code_start:code_end].strip() and not continuation:
          yield (code, comment, line_number )
          code       = ''
          comment    = ''
          continuing = False

    # Scan file for dependencies.
    #
    self._logger.logEvent( '  Scanning ' + sourceFilename )
    self._database.removeFile( sourceFilename )

    program_unit = None
    modules = []
    dependencies = []
    scope_stack = []
    pFUnitDriver = False
    for code, comment, line_number in lines_of_code( processed_source ):
      match = self._programPattern.match( code )
      if match:
        program_unit = match.group(1).lower()
        self._logger.logEvent( '    Contains program: ' + program_unit )
        self._database.addProgram( program_unit, sourceFilename )
        scope_stack.append( ('program', program_unit) )
        continue

      match = self._modulePattern.match( code )
      if match:
        program_unit = match.group( 1 ).lower()
        self._logger.logEvent( '    Contains module ' + program_unit )
        modules.append( program_unit )
        self._database.addModule( program_unit, sourceFilename )
        scope_stack.append( ('module', program_unit) )
        continue

      match = self._submodulePattern.match( code )
      if match:
        ancestorUnit = match.group(1)
        if ancestorUnit:
          ancestorUnit = ancestorUnit.lower()
        parentUnit  = match.group(2).lower()
        program_unit = match.group(3).lower()

        message = '{}Contains submodule {} of {}'.format( ' ' * 4,     \
                                                          program_unit, \
                                                          parentUnit )
        if ancestorUnit:
          message = message + '({})'.format( ancestorUnit )
        self._logger.logEvent( message )

        self._database.addModule( program_unit, sourceFilename )
        # I don't think it's necessary to append this to "modules".
        # It's really just a dependency.
        add_dependency( program_unit, parentUnit, True )
        scope_stack.append( ('submodule', program_unit) )
        continue

      match = self._subroutinePattern.match( code )
      if match and len(scope_stack) == 0: 
        # Only if this subroutine is a program unit.
        is_module    = match.group( 1 )
        program_unit = match.group( 2 ).lower()
        self._logger.logEvent( '    Contains subroutine ' + program_unit )
        modules.append( program_unit )
        self._database.add_procedure( program_unit, sourceFilename )
        scope_stack.append( ('subroutine', program_unit) )
        continue

      match = self._functionPattern.match( code )
      if match and len(scope_stack) == 0:
        # Only if this function is a program unit.
        is_module    = match.group( 1 )
        program_unit = match.group( 2 ).lower()
        self._logger.logEvent( '    Contains function ' + program_unit )
        modules.append( program_unit )
        self._database.add_procedure( program_unit, sourceFilename )
        scope_stack.append( ('function', program_unit) )
        continue

      match = self._endPattern.match( code )
      if match:
        end_scope = match.group(1)
        end_unit  = match.group(2)

        try:
          begin_scope, begin_unit = scope_stack[-1]

          if end_scope == begin_scope and end_unit == begin_unit:
            scope_stack.pop()
        except IndexError:
          message = 'Mismatched begin/end. Found "{end}" but stack empty'
          raise Exception( message.format( end=end_scope) )

      match = self._usePattern.match( code )
      if match is not None:
        moduleName = match.group( 1 ).lower()
        self._logger.logEvent( '    Depends on module ' + moduleName )
        add_dependency( program_unit, moduleName )
        continue

      match = self._externalPattern.match( code )
      if match is not None:
        moduleNames = [moduleName.strip()
                       for moduleName in match.group(1).split(',')]
        for moduleName in moduleNames:
          if moduleName is not None:
            self._logger.logEvent( '    Depends on external ' + moduleName )
            add_dependency( program_unit, moduleName )
        continue

      match = self._pFUnitPattern.match( code )
      if not pFUnitDriver and match is not None:
        self._logger.logEvent( '    Is driver' )
        pFUnitDriver = True

        includeFilename = os.path.join( os.path.dirname( sourceFilename ),
                                        'testSuites.inc' )
        with open( includeFilename, 'r' ) as includeFile:
          for line in includeFile:
            match = self._suitePattern.match( line )
            if match is not None:
              testGeneratorFunction = match.group( 1 )
              testModule = testGeneratorFunction.replace( '_suite', '' )
              self._logger.logEvent( '      Depends on module ' + testModule )
              self._database.addModuleCompileDependency( program_unit, testModule )
              self._database.addModuleLinkDependency( program_unit, testModule )
        continue

      for match in self._dependsPattern.finditer( comment ):
        name = match.group( 1 ).lower()
        extension = match.group( 2 )
        if name is not None:
            self._logger.logEvent( '    %s depends on call to %s ' % (program_unit, name) )
        add_dependency( program_unit, name )
        continue
