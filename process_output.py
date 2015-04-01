# imports:
import process_output_functions as fn
import numpy as np
import scipy as sci
import scipy.linalg as lin
import  matplotlib.cm as cm
import os

# gather inputs for creating "fake" files:
nx = int(raw_input("nx? \t\t\t"))
ny = int(raw_input("ny? \t\t\t"))
nz = int(raw_input("nz? \t\t\t"))
numParticles = int(raw_input("how many particles? \t"))
numFiles = int(raw_input("how many files? \t"))
simID = raw_input("simID? \t")
dim = -1
# nx = 100
# ny = 100
# nz = 100
# numParticles = 1000
# numFiles = 50

print("\n\nCreating fake files...")
# generate fake output files to test visualization functionality with:
if nz > 0:
    dim = 3
    fn.createTestOutputFiles_ALL3D(nx, ny, nz, numParticles, numFiles)
else:
    dim = 2
    fn.createTestOutputFiles_ALL2D(nx, ny, numParticles, numFiles)
print("*****************************\n")

# suck in all the gridpoint files and parse the text for the relevant data:
if(dim == 2):
    timesteps, E_min, E_max, B_min, B_max = fn.getGridData2D(numFiles)
elif(dim == 3):
    timesteps, E_min, E_max, B_min, B_max = fn.getGridData3D(numFiles)

# suck in all the particles files and parse the text for the relevant data:
if(dim == 2):
    timesteps, p_min, p_max, = fn.getParticleData2D(numFiles, timesteps)
elif (dim == 3):
    timesteps, p_min, p_max, = fn.getParticleData3D(numFiles, timesteps)

# normalize data:
timesteps = fn.normalizeData(timesteps, E_min, E_max, B_min, B_max, p_min, p_max)

# plot data:
cur_dir = os.getcwd()
path_name = "sim_" + str(simID)
path = os.path.join(cur_dir, path_name)

if(dim == 2):
    for t in xrange(0, numFiles):
        fn.plotDataForSingleTimeStep2D(timesteps[t]['gridpoints'], timesteps[t]['particles'], t, nx, ny, path)
elif(dim == 3):
    for t in xrange(0, numFiles):
        fn.plotDataForSingleTimeStep3D(timesteps[t]['gridpoints'], timesteps[t]['particles'], t, nx, ny, nz, path)

   
            