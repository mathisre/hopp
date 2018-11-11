from __future__ import print_function
import os
import sys
import numpy as np


try:
	T1 = float(sys.argv[1])
	T2 = float(sys.argv[2])
	dT = float(sys.argv[3])
except IndexError, ValueError:
	print("Put T1 , T2 and dT")
	exit()

Np = int(round((T2-T1)/dT))
print(Np)

T_lin = np.linspace(T1, T2, Np+1)
assert ((T_lin[1]-T_lin[0]) <= dT + 10**-4) and ((T_lin[1]-T_lin[0]) >= dT - 10**-4)
print("Temperatures:")
[print("{:.4f}, ".format(T_lin[i]), end='') for i in range(len(T_lin)-1)]
print("{:.4f}, ".format(T_lin[(len(T_lin)-1)]))



L = 100
n = 10
fill = 0.5
directory = "../dataMottNoTri/"
#directory = "quickData/"

for i in range(0,n):
	for T in T_lin:
		Ex = T/10
		Nt = int(1*10**7)
		Hz = 0
		seed = np.random.randint(10**7)
		seedconfig = np.random.randint(10**7)
		seedstate = np.random.randint(10**7)
		statement = "./glatzprogram outpre=%s%d_ Nt=%d Ex=%.7f Hz=%.7f T=%.7f seed2=%d rsconfig=%d rsstate=%d screen=0.5 fill=%.4f L=%d triJumps=0" % (directory,i,Nt,Ex,Hz,T,seed,seedconfig,seedstate,fill,L)
		print(statement)
		os.system(statement)
