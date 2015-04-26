#!/usr/bin/env python

# imports:
import process_output_functions as fn
import numpy as np
import scipy as sci
import scipy.linalg as lin
import matplotlib.cm as cm
import os


# gather inputs for creating "fake" files:
dim = int(raw_input("dim? \t"))
vidLen = int(raw_input("video length in secs? \t"))
simID = raw_input("simID? \t\t\t")
ext = "data"
win_or_linux = "win"

# grab params info from params file:
params_file = open('params.data', 'r')
if(dim == 2):
    params_file_lines = params_file.readlines() # params_file_lines is a list of the lines in params_file
    line = params_file_lines[1]
    nx = float((((line.split(","))[0]).split("="))[1])
    ny = float((((line.split(","))[1]).split("="))[1])
    line = params_file_lines[2]
    dx = float((((line.split(","))[0]).split("="))[1])
    dy = float((((line.split(","))[1]).split("="))[1])
    line = params_file_lines[3]
    numFiles = int((line.split("="))[1])
    line = params_file_lines[4]
    numParticles = int((line.split("="))[1])
if(dim == 3):
    params_file_lines = params_file.readlines() # params_file_lines is a list of the lines in params_file
    line = params_file_lines[1]
    nx = float((((line.split(","))[0]).split("="))[1])
    ny = float((((line.split(","))[1]).split("="))[1])
    nz = float((((line.split(","))[2]).split("="))[1])
    line = params_file_lines[2]
    dx = float((((line.split(","))[0]).split("="))[1])
    dy = float((((line.split(","))[1]).split("="))[1])
    dz = float((((line.split(","))[2]).split("="))[1])
    line = params_file_lines[3]
    numFiles = int((line.split("="))[1])
#    line = params_file_lines[4]
#    numParticles = int((line.split("="))[1])
params_file.close()

cur_dir = os.getcwd()
path_name = "sim_" + str(simID)
path = os.path.join(cur_dir, path_name)

# suck in all the gridpoint files and parse the text for the relevant data:
if(dim == 2):
    timesteps, E_min, E_max, B_min, B_max = fn.getGridData2D(numFiles, ext)
elif(dim == 3):
    timesteps, E_min, E_max, B_min, B_max = fn.getGridData3D(numFiles, ext, nx, ny, nz, dx, dy, dz)
    
# suck in all the particles files and parse the text for the relevant data:
if(dim == 2):
    timesteps, p_min, p_max, = fn.getParticleData2D(numFiles, timesteps, ext)
elif (dim == 3):
    timesteps, p_min, p_max, = fn.getParticleData3D(numFiles, timesteps, ext)
    
# normalize data:
timesteps = fn.normalizeData(timesteps, E_min, E_max, B_min, B_max, p_min, p_max, nx, ny, nz)
print("E_min is: " + str(E_min))
print("E_max is: " + str(E_max))
print("B_min is: " + str(B_min))
print("B_max is: " + str(B_max))

# plot data:    
if(dim == 2):
    for t in xrange(0, numFiles):
        fn.plotDataForSingleTimeStep2D(timesteps[t]['gridpoints'], timesteps[t]['particles'], t, nx, ny, dx, dy, path, 'B')
    for t in xrange(0, numFiles):
        fn.plotDataForSingleTimeStep2D(timesteps[t]['gridpoints'], timesteps[t]['particles'], t, nx, ny, dx, dy, path, 'E')
    
elif(dim == 3):
    for t in xrange(0, numFiles):
        fn.plotDataForSingleTimeStep3D(timesteps[t]['gridpoints'], timesteps[t]['particles'], t, nx, ny, nz, dx, dy, dz, path, 'B')
    for t in xrange(0, numFiles):
        fn.plotDataForSingleTimeStep3D(timesteps[t]['gridpoints'], timesteps[t]['particles'], t, nx, ny, nz, dx, dy, dz, path, 'E')

# elif makevid == 'y':
# if os.access(path, os.F_OK):
#     os.chdir(path)   
fn.makeVideos(path, 'B', win_or_linux, numFiles, vidLen)
fn.makeVideos(path, 'E', win_or_linux, numFiles, vidLen)

   
            
