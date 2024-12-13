#!/usr/bin/python
import os
from math import*
import time
os.system("cd ~/frw")
os.system("g++ -O3 frw.cpp MathAux.cpp -o frw.exe")
#os.system("mkdir /home/leo/frw/.msdbackup")


outdirname = "new"

L = 10000   # beware L_MAX @frw.cpp = 1048576    (  L = int(N/rho) )
os.system("mkdir " + outdirname)

N_repeats = 1 # 16/N_data_curves not included
N_data_curves = 1

for i in range(int(1/N_data_curves)): #for 16 cores
    rho = 0.8
    dt = 1000
    nfr = 1000
    q = 0.8
    for j in range(N_data_curves):
        N = int(L*rho)
        tmax = dt*nfr
        if os.fork() == 0:
            #os.system("xclock&")
            dirname = "/home/leo/frw/" + outdirname + "/msd_q" + str(q)
            dirname += "_rho" + str(rho) + "_dt" + str(dt) + "_nfr" + str(nfr)
            dirname += "_L" + str(L)

            os.system("mkdir " + dirname + " " + dirname + "/script")
            os.system("cp msd_run.py frw.cpp " + dirname + "/script/")
            dirname += "/data"
            os.system("mkdir " + dirname)
            if os.path.isfile(dirname + "/inputs.txt") == False:
                input_list = [q, rho, dt, nfr]
                inputs = open(dirname + "/inputs.txt", "w")
                for value in input_list:
                    inputs.write(str(value) + " ")
                inputs.close()
            for k in range(N_repeats):
                trialscnt = 1
                while(os.path.isdir(dirname + "/" + str(trialscnt))):
                    trialscnt += 1                
                 
                os.system("mkdir " + dirname + "/" + str(trialscnt))
            
                com = "nice -20 "
                com += "frw.exe " + str(q) + " " + str(N) + " " + str(L) + " " + str(tmax) + " " + str(nfr)  + " " + dirname + "/" + str(trialscnt) + "/ traj.dat;"
                os.system(com)
                print(com)
                os.system("xterm -e ~/bin/frw/msdruncal " + dirname + "/" + str(trialscnt) + " " + str(dt) + ";")
                
            exit()
        q += 0.2
        q = round(q, 4)
        #dt *= 0.6
        #dt = int(dt)
        #rho = round(rho, 2)
        time.sleep(3.1415)
