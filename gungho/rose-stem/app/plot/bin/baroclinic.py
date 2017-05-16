#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2017 Met Office. All rights reserved.
# For further details please refer to the file LICENCE which you should have
# received as part of this distribution.
##############################################################################
'''
Python script to plot xz slices along the y=0 and xy slices on a specified level 

This version takes nodal format output files and
interpolates onto a regular grid.

Filename hardcoded.

Levels are determined from the data

This version stitches together a directory of files
and extracts all levels so it can work in the serial
case where there is one file or the parallel case where
there is a file for each processor.

'''

import numpy as np
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.interpolate import griddata
import math
import glob
import sys

# Make an empty list to hold the levels we find in the data



def process_file_list(filestem, x, y, z, levels):

  # get the list of files to stitch together
  dirlist = glob.glob(filestem)

  # If no files are found then don't try to process them
  if len(dirlist) < 1:
    print("No files found to plot")
  else:

    for f in dirlist:
      print "processing file ", f
      fo = open(f, "r")

      # Step through all lines in the file, split the lines
      # and where the level matches the specifed one, append 
      # data to appropriate list 
      for strline in fo:
         strsplit = strline.split()
         # check we got a valid data line
         if (len(strsplit) == 5):
            # get the level
            level = float(strsplit[3])
            # Is the level already in the levels list?
            if (level in levels):
               # If it is then append the data into the correct list
               x[levels.index(level)].append(float(strsplit[0]))
               y[levels.index(level)].append(float(strsplit[1]))
               z[levels.index(level)].append(float(strsplit[4]))
            else:
               # add the level to the levels list and append
               # corresponding empty lists to x, y and z lists
               levels.append(level)
               x.append([])
               y.append([])
               z.append([])
               # ...and then append the data
               x[levels.index(level)].append(float(strsplit[0]))
               y[levels.index(level)].append(float(strsplit[1]))
               z[levels.index(level)].append(float(strsplit[4]))

      fo.close()
       
def make_figure(plotpath, field, timestep, x, y, z, levels):

   # get min and max of x,y data for plot axes
  xmin =  min(x[0])
  xmax = max(x[0])
  ymin =  min(y[0])
  ymax = max(y[0])
  zmin = min(levels)*1000.0
  zmax = max(levels)*1000.0

  r2d = 180.0/np.pi;
  nx,ny,nz = 360,180,len(levels)

  #create 2D plot
  x2d = np.linspace(xmin, xmax, nx)
  z2d = np.linspace(zmin, zmax, nz)
  y2d = np.linspace(ymin, ymax, ny)
  zi = np.zeros([ny,nx,len(levels)])
  xi, yi = np.meshgrid(x2d, y2d) 
  for p in xrange(len(levels)):
    zi[:,:,p] = griddata((np.asarray(x[p]), np.asarray(y[p])), np.asarray(z[p]), (xi, yi), method='linear')

  cc = np.linspace(np.amin(z),np.amax(z),13)

  # Pressure plot
  plotlevel = 1
  if field == 'theta':
    cc = np.linspace(220,330,12)
  elif field == 'exner':
    rd = 287.05
    p0 = 100000.0
    kappa = rd/1005.0
    zi = 0.01*zi**(1.0/kappa) * p0
    cc = np.linspace(900,1000,11)
    plotlevel = 0

  fig = plt.figure(figsize=(10,5))
  dz = zi[:,:,plotlevel]
  cf = plt.contourf(xi *r2d, yi * r2d, dz, cc)
  plt.colorbar(cf,  cmap=cm.spectral)
  cl = plt.contour(xi * r2d, yi * r2d, dz, cc, linewidths=0.5, colors='k')
  plt.xlabel('Longitude')
  plt.ylabel('Latitude')
  out_file_name = plotpath + "/" "slice_xy_" + field + "_" + timestep +  ".png"
  plt.savefig(out_file_name , bbox_inches='tight')
 
if __name__ == "__main__":

  try:
    config, datapath, fields, timesteps, plotpath = sys.argv[1:6]
  except ValueError:
    print("Usage: {0} <config> <datapath> <fields_list> <timestep_list> <plotpath>".format(sys.argv[0]))
    exit(1)

  # Split out the list of fields
  field_list = fields.split(':')

  # Split out the list of timesteps
  ts_list = timesteps.split(':')


  for ts in ts_list:
    for field in field_list:
      x      = []
      y      = []
      z      = []
      levels = []
      filestem =  datapath + "/" + config + "_nodal_" + field + "_" + ts + "*"      
      process_file_list(filestem, x, y, z, levels)

      # Only try to plot if we found some files for this timestep
      if len(levels) > 0:
        make_figure(plotpath, field, ts, x, y, z, levels)


