#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2017 Met Office. All rights reserved.
# For further details please refer to the file LICENCE which you should have
# received as part of this distribution.
##############################################################################

'''
Python script to plot all levels from a Dynamo output file. Levels are determined
from the data as they are different for different fields.

This version takes nodal format output files and interpolates onto a regular
grid.

This version stitches together a directory of files and extracts all levels
so it can work in the serial case where there is one file or the parallel
case where there is a file for each processor.

This version is for plotting under suites and accepts command line args
for the field and timestep to plot. It also plots to file rather than to screen 
'''
import numpy as np

# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import glob
import sys
from scipy.interpolate import griddata

levels = [] ; x = [] ; y = [] ; w3field = []

def process_file_list(filestem):

  # get the list of files to stitch together
  dirlist = glob.glob(filestem)

  for f in dirlist:

    fo = open(f, "r")

    # Step through all lines in the file, split the lines
    # and where the level matches the specifed one, append 
    # data to appropriate list 
    for strline in fo:
      strsplit = strline.split()
      if (len(strsplit) == 5): # check we got a valid data line

         level = float(strsplit[3]) # get the level

         if (level in levels):
            # If it is then append the data into the correct list
            x[levels.index(level)].append(float(strsplit[0]))
            y[levels.index(level)].append(float(strsplit[1]))
            w3field[levels.index(level)].append(float(strsplit[4]))
         else:
            # add the level to the levels list and append
            # corresponding empty lists to x, y and field lists
            levels.append(level)
            x.append([])
            y.append([])
            w3field.append([])
            # ...and then append the data
            x[levels.index(level)].append(float(strsplit[0]))
            y[levels.index(level)].append(float(strsplit[1]))
            w3field[levels.index(level)].append(float(strsplit[4]))

      if (len(strsplit) == 7):
         # get the level
         level = float(strsplit[3])

         dimension_to_plot = 4     # E/W direction ?
#         dimension_to_plot = 5      # N/S direction? 
#         dimension_to_plot = 6     # vertical direction ?

         if (level in levels):
#            # If it is then append the data into the correct list
            x[levels.index(level)].append(float(strsplit[0]))
            y[levels.index(level)].append(float(strsplit[1]))
            w3field[levels.index(level)].append(float(strsplit[dimension_to_plot]))
         else:
#            # add the level to the levels list and append
#            # corresponding empty lists to x, y and z lists
            levels.append(level)
            x.append([])
            y.append([])
            w3field.append([])
#            # ...and then append the data
            x[levels.index(level)].append(float(strsplit[0]))
            y[levels.index(level)].append(float(strsplit[1]))
            w3field[levels.index(level)].append(float(strsplit[dimension_to_plot]))

    fo.close()


def make_figure(field, timestep):

  fig = plt.figure(figsize=(15,10))
  plt.plot()

#  for p in xrange(len(levels)):
  for p in xrange(1):

    xmin, xmax, ymin, ymax =  min(x[p]), max(x[p]), min(y[p]), max(y[p])

    # Size of regular grid
    ny, nx = 200, 200 # Just set these to match the underlying ugrid mesh size for the biperiodic.

    # Generate a regular grid to interpolate the data onto
    xi, yi = np.meshgrid(np.linspace(xmin, xmax, nx), np.linspace(ymin, ymax, ny))
    # Interpolate using delaunay triangulation 
    w3fieldi = griddata((x[p], y[p]), w3field[p], (xi, yi), method='linear')

    cf = plt.pcolormesh(xi, yi, w3fieldi,cmap=cm.hsv,vmin=-1.0,vmax=10.0) # use this for biperiodic grid
    plt.colorbar(cf,cmap=cm.spectral)

    out_file_name = plotpath + "/" "biperiodic-cosmic_" + field + "_" + timestep +  ".png"
    plt.savefig(out_file_name , bbox_inches='tight')

  print out_file_name


if __name__ == "__main__":

  try:
    datapath, fields, timesteps, plotpath = sys.argv[1:5]
  except ValueError:
    print("Usage: {0} <datapath> <field_names> <timestep_list> <plotpath>".format(sys.argv[0]))
    exit(1)

  # Split out the list of fields
  field_list = fields.split(':')

  # Split out the list of timesteps
  ts_list = timesteps.split(':')

  for field in field_list:

    for ts in ts_list:

      # clear the lists in between plots
      del levels[:] ;  del x[:] ; del y[:] ; del w3field[:]
      filestem =  datapath + "/diagDynamo_nodal_" + field + "_" + ts + "*"

      process_file_list(filestem)
      # Only try to plot if we found some files for this timestep
      make_figure(field, ts)

    del levels[:] ; del x[:] ; del y[:] ; del w3field[:]

