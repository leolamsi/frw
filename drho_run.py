#!/usr/bin/python
import os
from math import*
import time
os.system("cd ~/frw")
os.system("g++ -O3 frw.cpp MathAux.cpp -o frw.exe")
#os.system("mkdir /home/leo/frw/.msdbackup")


L = 10000   # beware L_MAX @frw.cpp = 1048576    (  L = int(N/rho) )

N_repeats = 1 # repeat same set of inputs (serial on one thread)
N_rho = 1 # number of different rhos
N_q = 1 # number of different qs

q = 0.4
dt = 2000000
rho_init = 0.1
nfr = 500

add_q = 0
round_q_sigfig = 2
add_dt = 0
add_rho = 0
round_rho_sigfig = 4


dname = "new_runs"
os.system("mkdir " + dname)

for i in range(N_q): 
    rho = rho_init    
    qdname = dname + "/" + str(q)
    os.system("mkdir " + qdname)
    for j in range(N_rho):
        N = int(L*rho)
        tmax = dt*nfr
        if os.fork() == 0:
            dirname = "/home/leo/frw/" + qdname + "/msd_q" + str(q)
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
        rho += add_rho
        rho = round(rho, round_rho_sigfig)
        dt += add_dt
        time.sleep(3)
    q += add_q
    q = round(q, round_q_sigfig)
    #time.sleep(500)
