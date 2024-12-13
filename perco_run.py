import os
from math import*
import time
os.system('cd ~/frw')
os.system('g++ -O3 -fstack-protector-all perco_frw.cpp MathAux.cpp -o perco_frw.exe')
os.system('g++ -O3 -fstack-protector-all cal_sd.cpp -o cal_sd.exe')
#os.system('mkdir /home/leo/frw/.databackup')
#os.system('/home/leo/bin/renamedata')
#os.system('mkdir /home/leo/frw/data')


L = 500
dt = 10
nfr = 100
tmax = dt*nfr

q = 0.572
qstep = 0.001

N = 2
N_data_curves = 8
N_repeats = 3

outdirname = "perco_data"
os.system("mkdir /home/leo/frw/" + outdirname)

for j in range(N_data_curves):
    if os.fork() == 0:
        dirname = "/home/leo/frw/" + outdirname + "/" + str(N)
        os.system("mkdir " + dirname)
        dirname += "/dt" + str(dt)
        os.system("mkdir " + dirname)
        dirname += "/q" + str(q)
        os.system("mkdir " + dirname + " " + dirname + "/script")
        os.system("cp perco_run.py perco_frw.cpp cal_sd.cpp " + dirname + "/script/")
        if os.path.isfile(dirname + "/inputs.txt") == False:
            input_list = [N, q, dt]
            inputs = open(dirname + "/inputs.txt", "w")
            for value in input_list:
                inputs.write(str(value) + " ")
            inputs.close()
        for k in range(N_repeats):
            trialscnt = 1;
            while(os.path.isdir(dirname + "/" + str(trialscnt))):
                trialscnt += 1
            os.system("mkdir " + dirname + "/" + str(trialscnt))
                  
            com = "nice -20 "
            com += "perco_frw.exe " + str(q) + " " + str(N) + " " + str(L) + " " + str(tmax) + " " + str(nfr) + " " + dirname + "/" + str(trialscnt) + "/ raw; "
            com += "cal_sd.exe " + str(L) + " " + str(N) + " " + str(dt) + " " + dirname + "/" + str(trialscnt) + "/raw " + dirname + "/" + str(trialscnt) +  "/msd;"
            os.system(com)
            print(com)
                
        exit()
    q += qstep
    q = round(q, 4)
    time.sleep(3.1415)

