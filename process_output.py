#!/usr/bin/env python

# imports:
import process_output_functions as fn
import numpy as np
import scipy as sci
import scipy.linalg as lin
import matplotlib.cm as cm
import os


# gather inputs for creating "fake" files:
fakes = raw_input("create fake files? \t")
if(fakes == "y"):
    numParticles = int(raw_input("how many particles? \t"))
dim = raw_input("dim? \t")
vidLen = int(raw_input("video length in secs? \t"))
simID = raw_input("simID? \t\t\t")
ext = "data"
win_or_linux = "win"

# grab params info from params file:
with open("params.data", 'r') as params_file:
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

# # gather inputs for creating "fake" files:
# fakes = raw_input("create fake files? \t")
# ext = raw_input("extension? \t\t")
# nx = int(raw_input("nx? \t\t\t"))
# ny = int(raw_input("ny? \t\t\t"))
# nz = int(raw_input("nz? \t\t\t"))
# numParticles = int(raw_input("how many particles? \t"))
# numFiles = int(raw_input("how many files? \t"))
# vidLen = int(raw_input("video length in secs? \t"))
# win_or_linux = raw_input("windows or linux? \t")
# simID = raw_input("simID? \t\t\t")
# # makevid = raw_input("make video? (y,n) \t")

# nx = 100
# ny = 100
# nz = 100
# numParticles = 1000
# numFiles = 50
# vidLen = 10
# win_or_linux = "win"

cur_dir = os.getcwd()
path_name = "sim_" + str(simID)
path = os.path.join(cur_dir, path_name)

if(fakes == 'y'):
    print("\n\nCreating fake files...")
    # generate fake output files to test visualization functionality with:
    if nz > 0:
        dim = 3
        fn.createTestOutputFiles_ALL3D(nx, ny, nz, numParticles, numFiles, ext)
    else:
        dim = 2
        fn.createTestOutputFiles_ALL2D(nx, ny, numParticles, numFiles, ext)
    print("*****************************\n")
    
# suck in all the gridpoint files and parse the text for the relevant data:
if(dim == 2):
    timesteps, E_min, E_max, B_min, B_max = fn.getGridData2D(numFiles, ext)
elif(dim == 3):
    timesteps, E_min, E_max, B_min, B_max = fn.getGridData3D(numFiles, ext)
    
# suck in all the particles files and parse the text for the relevant data:
if(dim == 2):
    timesteps, p_min, p_max, = fn.getParticleData2D(numFiles, timesteps, ext)
elif (dim == 3):
    timesteps, p_min, p_max, = fn.getParticleData3D(numFiles, timesteps, ext)
    
# normalize data:
timesteps = fn.normalizeData(timesteps, E_min, E_max, B_min, B_max, p_min, p_max)
    
# plot data:    
if(dim == 2):
    for t in xrange(0, numFiles):
        fn.plotDataForSingleTimeStep2D(timesteps[t]['gridpoints'], timesteps[t]['particles'], t, nx, ny, path, 'B')
    for t in xrange(0, numFiles):
        fn.plotDataForSingleTimeStep2D(timesteps[t]['gridpoints'], timesteps[t]['particles'], t, nx, ny, path, 'E')
    
elif(dim == 3):
    for t in xrange(0, numFiles):
        fn.plotDataForSingleTimeStep3D(timesteps[t]['gridpoints'], timesteps[t]['particles'], t, nx, ny, nz, path, 'B')
    for t in xrange(0, numFiles):
        fn.plotDataForSingleTimeStep3D(timesteps[t]['gridpoints'], timesteps[t]['particles'], t, nx, ny, nz, path, 'E')

# elif makevid == 'y':
# if os.access(path, os.F_OK):
#     os.chdir(path)   
fn.makeVideos(path, 'B', win_or_linux, numFiles, vidLen)
fn.makeVideos(path, 'E', win_or_linux, numFiles, vidLen)

   
            
