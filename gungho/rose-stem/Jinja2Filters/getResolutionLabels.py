#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Implements a Jinja2 filter to generate strings for resolution labels/options.
'''
from jinja2 import contextfilter

@contextfilter
def getResolutionLabels(contextfilter, resolution):
    '''
    Returns labels used by Jinja2 to generate app information for
    mesh resolution and timestep options

    @param [in] context    Jinja2 instance to run macro against.
    @param [in] resolution List in format where 1st entry is always spatial
                           the option code for the mesh resolution, followed by
                           the timestep resolutions to run at (if any)
    @return Mesh resolution option
            Time resolution value(s),
            Mesh resolution label,
            Time resolution labels,
            Rose mesh resolution option string
    '''
    
    mesh_resolution_label  = ''
    time_resolution_labels = []
    mesh_resolution_option = ''
    rose_resolution_option_str = ''

    # Set mesh values/labels
    if resolution[0] != 'default':
        mesh_resolution_option = resolution[0]

    if mesh_resolution_option != '':
        mesh_resolution_label = "_" + mesh_resolution_option
        rose_resolution_option_str = '--opt-conf-key=' + mesh_resolution_option

    # Set timestep values/labels
    if len(resolution) > 1:
      time_resolution_values = resolution[1:]
    else:
      time_resolution_values = ['']

    for timestep in time_resolution_values:
        if timestep != '':
            time_resolution_labels.append('_dt-' + str(timestep).replace('.','p'))
        else:
            time_resolution_labels.append('')

    return mesh_resolution_option, time_resolution_values, mesh_resolution_label, \
           time_resolution_labels, rose_resolution_option_str
