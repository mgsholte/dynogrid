# imports:
import random

def createTestOutputFiles():
    itNum = 0
    dt = 0.000000000000000018
    # create a fake output file for itNum_grid.txt
    f_name = str(itNum) + "_grid.txt"
    with open(f_name, 'a+') as test_file:
        test_file.write("GRIDPOINTS\n")
        test_file.write("itNum=" + str(itNum) + "\n")
        test_file.write("Time=" + str(itNum*dt) + "\n")
        test_file.write("TimeStep=" + str(dt) + "\n")
        test_file.write("GridSize:nx=100,ny=100\n")
        
        for i in xrange(0,100):
            for j in xrange(0,100):       
                if itNum == 0:
                    test_file.write("coordinates:"+str(i)+","+str(j)+","+"Ex=0,Ey=0,Ez=0,Bx=0,By=0,Bz=0\n")
                else:
                    Ex = random.random()/1000.0
                    Ey = random.random()/1000.0
                    Ez = random.random()/1000.0
                    Bx = random.random()/1000.0
                    By = random.random()/1000.0
                    Bz = random.random()/1000.0
                    test_file.write("coordinates:" + str(i) + "," + str(j) + "," + "Ex=" + str(Ex) + "," + "Ey=" + str(Ey) + "," + "Ez=" + str(Ez) + "," + "Bx=" + str(Bx) + "," + "By=" + str(By) + "," + "Bz=" + str(Bz) + "\n")
    
    # create a fake output file for itNum_particles.txt
    f_name2 = str(itNum) + "_particles.txt"
    with open(f_name2, 'a+') as test_file2:
        test_file2.write("PARTICLES\n")
        test_file2.write("itNum=" + str(itNum) + "\n")
        test_file2.write("Time=" + str(itNum*dt) + "\n")
        test_file2.write("TimeStep=" + str(dt) + "\n")
        test_file2.write("GridSize:nx=100,ny=100\n")
        
        for i in xrange(0,1000):
                if itNum == 0:
                    pos_x = (random.random()*10.0)+45.0
                    pos_y = (random.random()*10.0)+45.0
                    test_file2.write("pos_x=" + str(pos_x) + "," + "pos_y=" + str(pos_y) + "," + "p_x=0,p_y=0,p_z=0\n")
                else:
                    px = random.random()/1000000000.0
                    py = random.random()/1000000000.0
                    pz = random.random()/1000000000.0
                    pos_x = (random.random()*10.0)+45.0
                    pos_y = (random.random()*10.0)+45.0
                    test_file2.write("pos_x=" + str(pos_x) + "," + "pos_y=" + str(pos_y) + "," + "p_x=" + str(px) + "," + "p_y=" + str(py) + "," + "p_z=" + str(pz) + "\n")
        