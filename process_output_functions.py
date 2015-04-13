# imports:
import random
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def createTestOutputFile2D(itNum, nx, ny, numParticles, numFiles, ext):
    dt = 0.000000000000000018
    # create a fake output file for itNum_grid.txt
    f_name = str(itNum) + "_grid." + str(ext)
    with open(f_name, 'w+') as test_file:
        test_file.write("GRIDPOINTS\n")
        test_file.write("itNum=" + str(itNum) + "\n")
        test_file.write("Time=" + str(itNum*dt) + "\n")
        test_file.write("TimeStep=" + str(dt) + "\n")
        test_file.write("GridSize:nx=100,ny=100\n")
        
        for i in xrange(0, ny):
            for j in xrange(0, nx):       
                if itNum == 0:
                    test_file.write(str(j)+","+str(i)+","+ "0.0,0.0\n")
                else:
                    E = random.random()/1000.00
                    B = random.random()/1000.0
                    test_file.write(str(j) + "," + str(i) + "," + str(E) + "," + str(B) + "\n")
    
    # create a fake output file for itNum_particles.txt
    f_name2 = str(itNum) + "_particles." + str(ext)
    with open(f_name2, 'w+') as test_file2:
        test_file2.write("PARTICLES\n")
        test_file2.write("itNum=" + str(itNum) + "\n")
        test_file2.write("Time=" + str(itNum*dt) + "\n")
        test_file2.write("TimeStep=" + str(dt) + "\n")
        test_file2.write("GridSize:nx=100,ny=100\n")
        
        for i in xrange(0, numParticles):
                if itNum == 0:
                    pos_x = (random.random()*(nx/10.0))+((nx/2.0)-5.0)
                    pos_y = (random.random()*(ny/10.0))+((ny/2.0)-5.0)
                    test_file2.write("ptcl:" + str(i) + "," + str(pos_x) + "," + str(pos_y) + "," + "0.0\n")
                else:
                    p = random.random()/1000000000.0
                    multFac = (float(itNum)/float(numFiles))
                    if multFac < 0.1:
                        multFac = 0.1
                    xrng = nx*multFac
                    xspread = (nx-xrng)/2.0
                    yrng = ny*multFac
                    yspread = (ny-yrng)/2.0
                    pos_x = (random.uniform(xspread, nx-xspread))
                    pos_y = (random.uniform(yspread, ny-yspread))
                    test_file2.write("ptcl:" + str(i) + "," + str(pos_x) + "," + str(pos_y) + "," + str(p) + "\n")

def createTestOutputFile3D(itNum, nx, ny, nz, numParticles, numFiles, ext):
    dt = 0.000000000000000018
    # create a fake output file for itNum_grid.txt
    f_name = str(itNum) + "_grid." + str(ext)
    with open(f_name, 'w+') as test_file:
        test_file.write("GRIDPOINTS\n")
        test_file.write("itNum=" + str(itNum) + "\n")
        test_file.write("Time=" + str(itNum*dt) + "\n")
        test_file.write("TimeStep=" + str(dt) + "\n")
        test_file.write("GridSize:nx=100,ny=100,nz=100\n")
        
        for i in xrange(0, ny):
            for j in xrange(0, nx): 
                for k in xrange(0, nz):      
                    if itNum == 0:
                        test_file.write(str(j)+","+ str(i) +"," + str(k) + ",0.0,0.0\n")
                    else:
                        E = random.random()/1000.0
                        B = random.random()/1000.0
                        test_file.write(str(j) + "," + str(i) + "," + str(k) + "," + str(E) + "," + str(B) + "\n")
    
    # create a fake output file for itNum_particles.txt
    f_name2 = str(itNum) + "_particles." + str(ext)
    with open(f_name2, 'w+') as test_file2:
        test_file2.write("PARTICLES\n")
        test_file2.write("itNum=" + str(itNum) + "\n")
        test_file2.write("Time=" + str(itNum*dt) + "\n")
        test_file2.write("TimeStep=" + str(dt) + "\n")
        test_file2.write("GridSize:nx=100,ny=100\n")
        
        for i in xrange(0, numParticles):
                if itNum == 0:
                    pos_x = (random.random()*(nx/10.0))+((nx/2.0)-5.0)
                    pos_y = (random.random()*(ny/10.0))+((ny/2.0)-5.0)
                    pos_z = (random.random()*(nz/10.0))+((nz/2.0)-5.0)
                    test_file2.write("ptcl:" + str(i) + "," + str(pos_x) + "," + str(pos_y) + "," + str(pos_z) + ",0.0\n")
                else:
                    p = random.random()/1000000000.0
                    multFac = (float(itNum)/float(numFiles))
                    if multFac < 0.1:
                        multFac = 0.1
                    xrng = nx*multFac
                    xspread = (nx-xrng)/2.0
                    yrng = ny*multFac
                    yspread = (ny-yrng)/2.0
                    zrng = nz*multFac
                    zspread = (nz-zrng)/2.0
                    pos_x = (random.uniform(xspread, nx-xspread))
                    pos_y = (random.uniform(yspread, ny-yspread))
                    pos_z = (random.uniform(zspread, nz-zspread))
                    test_file2.write("ptcl:" + str(i) + "," + str(pos_x) + "," + str(pos_y) + "," + str(pos_z) + "," + str(p) + "\n")


def createTestOutputFiles_ALL2D(nx, ny, numParticles, numFiles, ext):
    for i in xrange(0, numFiles):
        createTestOutputFile2D(i, nx, ny, numParticles, numFiles, ext)
        
def createTestOutputFiles_ALL3D(nx, ny, nz, numParticles, numFiles, ext):
    for i in xrange(0, numFiles):
        createTestOutputFile3D(i, nx, ny, nz, numParticles, numFiles, ext)        

def getGridData2D(numFiles, ext):   
    timesteps = list()
    # suck in all the gridpoint files and parse the text for the relevant data:
    E_min = 0.0
    E_max= 0.0
    B_min = 0.0
    B_max = 0.0
    fname_grid = "_grid." + str(ext)
    for i in xrange(0, numFiles):
        fname = str(i) + str(fname_grid)
        gridpoints = dict()
        Xs = list()
        Ys = list()
        Es = list()
        Bs = list()
        with open(fname, 'r') as grid_file:
            grid_file_lines = grid_file.readlines() # grid_file_lines is a list of the lines in grid_file
            line_ct = 0
            for line in grid_file_lines:
                line_ct += 1
                if line_ct <= 5:
                    continue
                # parse eachline to get needed data:
                eachline = line.split(',') # eachline is a list of the results of the split
                xcoord = float(eachline[0])
                ycoord = float(eachline[1])
                E = float(eachline[2])
                B = float(eachline[3])
                if E < E_min:
                    E_min = E
                if E > E_max:
                    E_max = E
                if B < B_min:
                    B_min = B
                if B > B_max:
                    B_max = B
                Xs.append(xcoord)
                Ys.append(ycoord)
                Es.append(E)
                Bs.append(B)
        # store data for gridpoints in a dictionary:
        gridpoints = {'Xs':Xs, 'Ys':Ys, 'Es':Es, 'Bs':Bs}
        timesteps.append({'gridpoints':gridpoints,'particles':''})
    return timesteps, E_min, E_max, B_min, B_max

def getParticleData2D(numFiles, timesteps, ext):
    p_min = 0.0
    p_max = 0.0        
    # suck in all the particles files and parse the text for the relevant data:            
    fname_particles = "_particles." + str(ext)
    for i in xrange(0, numFiles):
        fname = str(i) + str(fname_particles)
        particles = dict()
        Xs = list()
        Ys = list()
        ps = list()
        with open(fname, 'r') as particle_file:
            particles = list()
            particle_file_lines = particle_file.readlines() # grid_file_lines is a list of the lines in grid_file
            line_ct = 0
            for line in particle_file_lines:
                line_ct += 1
                if line_ct <= 5:
                    continue
                # parse eachline to get needed data:
                eachline = line.split(',') # eachline is a list of the results of the split
                xcoord = float(eachline[1])
                ycoord = float(eachline[2])
                p = float(eachline[3])
                if p < p_min:
                    p_min = p
                if p > p_max:
                    p_max = p
                # store data for in a dictionary:
                Xs.append(xcoord)
                Ys.append(ycoord)
                ps.append(p)
        particles = {'Xs':Xs, 'Ys':Ys, 'ps':ps}
        (timesteps[i])['particles'] = particles
    return timesteps, p_min, p_max


def getGridData3D(numFiles, ext, nx, ny, nz):
    numLines = int(nx*ny*nz)
    timesteps = list()
    # suck in all the gridpoint files and parse the text for the relevant data:
    E_min = 0.0
    E_max= 0.0
    B_min = 0.0
    B_max = 0.0
    fname_grid = "_grid." + str(ext)
    for i in xrange(0, numFiles):
        fname = str(i) + str(fname_grid)
        
        gridpoints = dict()
        # Xs = list()
        # Ys = list()
        # Zs = list()
        # Bs = list()
        # Es = list()
        Xs = np.empty(numLines, dtype=np.float64)
        Ys = np.empty(numLines, dtype=np.float64)
        Zs = np.empty(numLines, dtype=np.float64)
        Bs = np.empty(numLines, dtype=np.float64)
        Es = np.empty(numLines, dtype=np.float64)
        with open(fname, 'r') as grid_file:
            # grid_file_lines = grid_file.readlines() # grid_file_lines is a list of the lines in grid_file
            line_ct = -1
            for line in grid_file:
                line_ct += 1
                if line_ct < 5:
                    continue
                # parse eachline to get needed data:
                eachline = line.split(',') # eachline is a list of the results of the split
                xcoord = float(eachline[0])
                ycoord = float(eachline[1])
                zcoord = float(eachline[2])
                E = float(eachline[3])
                B = float(eachline[4])
                # print("E is: " + str(E) + "\tB is: " + str(B))
                if E < E_min:
                    E_min = E
                if E > E_max:
                    E_max = E
                if B < B_min:
                    B_min = B
                if B > B_max:
                    B_max = B
                # Xs.append(xcoord)
                # Ys.append(ycoord)
                # Zs.append(zcoord)
                # Es.append(E)
                # Bs.append(B)
                Xs[line_ct - 5] = xcoord
                Ys[line_ct - 5] = ycoord
                Zs[line_ct - 5] = zcoord
                Es[line_ct - 5] = E
                Bs[line_ct - 5] = B
                # store data in a dictionary:
        gridpoints = {'Xs':Xs, 'Ys':Ys, 'Zs':Zs, 'Es':Es, 'Bs':Bs}
        timesteps.append({'gridpoints':gridpoints,'particles':''})
    return timesteps, E_min, E_max, B_min, B_max
        
def getParticleData3D(numFiles, timesteps, ext, numParticles):
    p_min = 0.0
    p_max = 0.0        
    # suck in all the particles files and parse the text for the relevant data:            
    fname_particles = "_particles." + str(ext)
    for i in xrange(0, numFiles):
        fname = str(i) + str(fname_particles)
        particles = dict()
        Xs = np.empty(numParticles, dtype=np.float64)
        Ys = np.empty(numParticles, dtype=np.float64)
        Zs = np.empty(numParticles, dtype=np.float64)
        ps = np.empty(numParticles, dtype=np.float64)
        with open(fname, 'r') as particle_file:
            particles = list()
            particle_file_lines = particle_file.readlines() # grid_file_lines is a list of the lines in grid_file
            line_ct = -1
            for line in particle_file_lines:
                line_ct += 1
                if line_ct < 5:
                    continue
                # parse eachline to get needed data:
                eachline = line.split(',') # eachline is a list of the results of the split
                xcoord = float(eachline[1])
                ycoord = float(eachline[2])
                zcoord = float(eachline[3])
                p = float(eachline[4])
                # print("p_xcoord is: " + str(xcoord) + "\tp_ycoord is: " + str(ycoord) + "\tp_zcoord is: " + str(zcoord))
                # print("p is: " + str(p))
                if p < p_min:
                    p_min = p
                if p > p_max:
                    p_max = p
                # store data in a dictionary:
                Xs[line_ct-5] = xcoord
                Ys[line_ct-5] = ycoord
                Zs[line_ct-5] = zcoord
                ps[line_ct-5] = p
        particles = {'Xs':Xs, 'Ys':Ys, 'Zs':Zs, 'ps':ps}
        (timesteps[i])['particles'] = particles
    return timesteps, p_min, p_max
        
        
def normalizeData(timesteps, E_min, E_max, B_min, B_max, p_min, p_max, nx, ny, nz, numParticles):
    # (x-min)/(max-min)
    for t in xrange(0, len(timesteps)):
        # normalize gridpoint data:
        # for grdpt in xrange(0, len(((timesteps[t])['gridpoints'])['Es'])):
        for grdpt in xrange(0, int(nx*ny*nz)):
            # normalize E and B in each gridpoint dictionary:
            if E_min == E_max:
                (((timesteps[t])['gridpoints'])['Es'])[grdpt] = 0.5
                # print("NOT GOOD! E_min == E_max")
            else:
                E_temp = (((timesteps[t])['gridpoints'])['Es'])[grdpt]
                (((timesteps[t])['gridpoints'])['Es'])[grdpt] = ((E_temp-E_min)/(E_max-E_min))
                # print("E_normed is: " + str((((timesteps[t])['gridpoints'])['Es'])[grdpt]))
            if B_min == B_max:
                (((timesteps[t])['gridpoints'])['Bs'])[grdpt] = 0.5
                # print("NOT GOOD! B_min == B_max")
            else:
                B_temp = (((timesteps[t])['gridpoints'])['Bs'])[grdpt]
                (((timesteps[t])['gridpoints'])['Bs'])[grdpt] = ((B_temp-B_min)/(B_max-B_min))
                # print("B_normed is: " + str((((timesteps[t])['gridpoints'])['Bs'])[grdpt]))
        # normalize particle data:        
        # for ptcl in xrange(0, len(((timesteps[t])['particles'])['ps'])):
        for ptcl in xrange(0, numParticles):
            # normalize p in each particle dictionary:
            if p_min == p_max:
                (((timesteps[t])['particles'])['ps'])[ptcl] = 0.5
                # print("p_normed is: " + str((((timesteps[t])['particles'])['ps'])[ptcl]))
            else:
                p_temp = (((timesteps[t])['particles'])['ps'])[ptcl]
                (((timesteps[t])['particles'])['ps'])[ptcl] = ((p_temp-p_min)/(p_max-p_min))
                # print("p_normed is: " + str((((timesteps[t])['particles'])['ps'])[ptcl]))
    return timesteps



def plotDataForSingleTimeStep2D(gridpoints, particles, itNum, nx, ny, dx, dy, path, E_or_B):
    itNum = '{0:05d}'.format(itNum)
    if E_or_B == 'B':
        title = "B_field_" + str(itNum)
    elif E_or_B == 'E':
        title = "E_field_" + str(itNum)

    fig = plt.figure(figsize=(12,7))
    ax1 = fig.add_subplot(111)
    gridpoints_Xs = np.array(gridpoints['Xs'])
    gridpoints_Ys = np.array(gridpoints['Ys'])
    gridpoints_Es = np.array(gridpoints['Es'])
    gridpoints_Bs = np.array(gridpoints['Bs'])

    fig.suptitle(title, fontsize='20')
    SM = cm.ScalarMappable()
    cMap = SM.get_cmap()
    plt.xlabel('X', fontsize=16)
    plt.ylabel('Y', fontsize=16)
    plt.xlim(0, nx)
    plt.ylim(0, ny)
    if E_or_B == 'B':
        cb1 = ax1.scatter(gridpoints_Xs, gridpoints_Ys, c=gridpoints_Bs, marker=u'o', cmap='binary', linewidths=0, alpha=.1, label='B Field')
    elif E_or_B == 'E':
        cb1 = ax1.scatter(gridpoints_Xs, gridpoints_Ys, c=gridpoints_Es, marker=u'o', cmap='binary', linewidths=0, alpha=.3, label='E Field')

    particles_Xs = np.array(particles['Xs'])
    particles_Ys = np.array(particles['Ys'])
    particles_ps = np.array(particles['ps'])
    cb2 = ax1.scatter(particles_Xs, particles_Ys, c=particles_ps, marker=u'o', cmap='coolwarm', linewidths=.3, label='Particles, mapped by momentum')

    legend = ax1.legend(loc='upper right', shadow=True)
    plt.colorbar(cb1)
    plt.colorbar(cb2)
#     plt.show()
    if os.access(path, os.F_OK):
        os.chdir(path)
    else:
        os.mkdir(path)
        os.chdir(path)
    fig.savefig(title)
    plt.close(fig)
    
    
def plotDataForSingleTimeStep3D(gridpoints, particles, itNum, nx, ny, nz, dx, dy, dz, path, E_or_B):
    itNum = '{0:05d}'.format(itNum)
    if E_or_B == 'B':
        title = "B_field_" + str(itNum)
    elif E_or_B == 'E':
        title = "E_field_" + str(itNum)    

    fig = plt.figure(figsize=(12,7))
    ax1 = fig.add_subplot(111, projection='3d')
    # gridpoints_Xs = np.array(gridpoints['Xs'])
    # gridpoints_Ys = np.array(gridpoints['Ys'])
    # gridpoints_Zs = np.array(gridpoints['Zs'])
    # gridpoints_Es = np.array(gridpoints['Es'])
    # gridpoints_Bs = np.array(gridpoints['Bs'])
    gridpoints_Xs = gridpoints['Xs']
    gridpoints_Ys = gridpoints['Ys'] # HACK (not anymore)
    gridpoints_Zs = gridpoints['Zs'] # HACK (not anymore)
    gridpoints_Es = gridpoints['Es']
    gridpoints_Bs = gridpoints['Bs']

    fig.suptitle(title, fontsize='20')
    # SM = cm.ScalarMappable()
    # cMap = SM.get_cmap()
    myCmap = cm.get_cmap("binary")
    myCmap.set_under('w', 0)

                        
    ax1.set_xlabel('X', fontsize=16)
    ax1.set_ylabel('Y', fontsize=16)
    ax1.set_zlabel('Z', fontsize=16)
    ax1.set_xlim(0, nx*dx)
    ax1.set_ylim(0, ny*dy)
    ax1.set_zlim(0, nz*dz)
    
    if E_or_B == 'B':
        cb1 = ax1.scatter(gridpoints_Xs, gridpoints_Ys, gridpoints_Zs, c=gridpoints_Bs, marker=u'.', cmap=myCmap, vmin=.2, linewidths=0, alpha=.1, label='B Field')
    elif E_or_B == 'E':
        cb1 = ax1.scatter(gridpoints_Xs, gridpoints_Ys, gridpoints_Zs, c=gridpoints_Es, marker=u'.', cmap=myCmap, vmin=.5, linewidths=0, alpha=.1, label='E Field')

    # particles_Xs = particles['Xs']
    # particles_Ys = particles['Ys']
    # particles_Zs = particles['Zs']
    # particles_ps = particles['ps']
    # for i in xrange(0, len(particles['Xs'])):
        # print("particles['Xs']["+str(i)+"] is: " + str(particles['Xs'][i]))
        # print("particles['Ys']["+str(i)+"] is: " + str(particles['Ys'][i]))
        # print("particles['Zs']["+str(i)+"] is: " + str(particles['Zs'][i]))
    cb2 = ax1.scatter(particles['Xs'], particles['Ys'], particles['Zs'], c=particles['ps'], marker=u'o', cmap='coolwarm', linewidths=.3, label='Particles, mapped by momentum')
   
    legend = ax1.legend(loc='upper right', shadow=True)
    plt.colorbar(cb1)
    plt.colorbar(cb2)

    #plt.show()
    if os.access(path, os.F_OK):
        os.chdir(path)
    else:
        os.mkdir(path)
        os.chdir(path)
    fig.savefig(title)
    plt.close()

    
def makeVideos(path, E_or_B, win_or_linux, numFiles, vidLen):
    img_names = ""
    fname = ""
    fps = int(float(numFiles)/float(vidLen))
    
    if win_or_linux == "win":
        if E_or_B == 'B':
            img_names = "B_field_%5d.png"
            fname = "B_field.mp4"
            command = "ffmpeg -f image2 -r " + str(fps) + " -i " + str(img_names) + " -vcodec mpeg4 -y " + str(fname)
            os.system(command)
        elif E_or_B == 'E':
            img_names = "E_field_%5d.png"
            fname = "E_field.mp4"
            command = "ffmpeg -f image2 -r " + str(fps) + " -i " + str(img_names) + " -vcodec mpeg4 -y " + str(fname)
            os.system(command)
    elif win_or_linux == "linux":
        if E_or_B == 'B':
            img_names = "B_field_%5d.png"
            fname = "B_field.mp4"
#             command = "ffmpeg -f image2 -r " + str(fps) + " -i " + str(img_names) + " -vcodec mpeg4 -y " + str(fname)
#             os.system(command)
        elif E_or_B == 'E':
            img_names = "E_field_%5d.png"
            fname = "E_field.mp4"
#             command = "ffmpeg -f image2 -r " + str(fps) + " -i " + str(img_names) + " -vcodec mpeg4 -y " + str(fname)
#             os.system(command)
    
    
