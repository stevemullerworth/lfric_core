#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a Jinja2 filter to strip the resolutions info.
'''
from jinja2 import contextfilter
import re
import ast
import science_parser

@contextfilter
def getResolutionMacro(context, call):
    '''
    Takes a string return list of resolutions to use
        crun: Number of runs to do in the crun.
        ...
    @param [in] context Jinja2 instance to run macro against.
    @return List resulting from retrieving the resolutions.
    '''

    if call.find('(') == -1:
        macroName = call
        arguments = []
    else:
        macroName = call[:call.index('(')]
        arguments = science_parser.science_parser(call[call.index('(')+1:call.rindex(')')])

    normalArguments  = [argument for argument in arguments \
                        if argument.find('=') == -1]
    keywordArguments = [argument for argument in arguments \
                        if argument.find('=') != -1]

    argumentList = []
    for argument in normalArguments:
        # Remove quote marks (if present) from the
        # argument string
        if argument[0] == '"':
            arg=argumentList.append( argument[1:-1] )
        else:
            arg=argumentList.append( argument )

    if len(argumentList) >= 2:

        app_name      = argumentList[0] # First argument is always the app_name
        configuration = argumentList[1] # Second argument is always the configuration
        appKey = app_name + '_' + configuration

        argumentDictionary = {}
        for argument in keywordArguments:
            key, value = re.split(' *= *', argument)
            argumentDictionary[key] = value

        # Return info about the resolutions arguments
        resList=[]
        if 'resolutions' in argumentDictionary.keys():
            resList = ast.literal_eval(argumentDictionary['resolutions'])

        resSupportMeshes=['']
        if 'support_meshes' in argumentDictionary.keys():
            resSupportMeshes = ast.literal_eval(argumentDictionary['support_meshes'])

        # We now consider the possibility that entries for the resolutions
        # can either just be the resolution (i.e. a string) or a
        # (resolution,timestep) pair.  For the latter we generate a
        # dictionary relating the resolution to the timestep.
        resDict={}
        for entry in resList:
            if type(entry) in [type(()),type([])]:
                if entry[0] not in resDict.keys():
                    resDict[entry[0]]=set()
                if len(entry) > 1:
                    if type(entry[1]) in [type(1),type(1.0)]:
                        resDict[entry[0]].add(entry[1])
                    if type(entry[1]) in [type(()),type([])]:
                        resDict[entry[0]].update(entry[1])

        return_value = appKey, resList, resDict, resSupportMeshes

    else:

        return_value = None, None, None, None

    return return_value
