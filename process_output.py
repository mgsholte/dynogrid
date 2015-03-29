# imports:
import random

# create a fake output file for itNum_grid.txt
with open("0_grid.txt", 'a+') as test_file:
    for i in xrange(0,100):
        for j in xrange(0,100):
            if i == 0 and j == 0:
                test_file.write("coordinates:"+str(i)+","+str(j)+","+"Ex=0,Ey=0,Ez=0,Bx=0,By=0,Bz=0\n")
            else:
                Ex = random.random()/1000.0
                Ey = random.random()/1000.0
                Ez = random.random()/1000.0
                Bx = random.random()/1000.0
                By = random.random()/1000.0
                Bz = random.random()/1000.0
                test_file.write("coordinates:" + str(i) + "," + str(j) + "," + "Ex=" + str(Ex) + "," + "Ey=" + str(Ey) + "," + "Ez=" + str(Ez) + "," + "Bx=" + str(Bx) + "," + "By=" + str(By) + "," + "Bz=" + str(Bz))

# create a fake output file for itNum_particles.txt
with open("0_particles.txt", 'a+') as test2_file:
    for i in xrange(0,1000):
            if i == 0:
                pos_x = (random.random()*10.0)+45.0
                pos_y = (random.random()*10.0)+45.0
                test2_file.write("pos_x=" + str(pos_x) + "," + "pos_y=" + str(pos_y) + "," + "p_x=0,p_y=0,p_z=0")
            else:
                px = random.random()/1000000000.0
                py = random.random()/1000000000.0
                pz = random.random()/1000000000.0
                pos_x = (random.random()*10.0)+45.0
                pos_y = (random.random()*10.0)+45.0
                test2_file.write("pos_x=" + str(pos_x) + "," + "pos_y=" + str(pos_y) + "," + "p_x=" + str(px) + "," + "p_y=" + str(py) + "," + "p_z=" + str(pz))
                








