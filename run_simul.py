import os
import sys
import numpy as np

T_lin = np.linspace(0.01, 1, 21)
def timesteps(T):
	return (10000/T**0.5)

E_lin = np.linspace(0, 1, 21)

for T_ in T_lin:
	for Ex in E_lin:
		Nt = int(round(timesteps(T_)))
		T = round(T_*100)/100.
		print("./glatzprogram triJumps=0 outpre=../../data/noTri/ writelines=500 Nt=" + str(Nt) + " Ex=" + str(Ex) + " Hz=0 " + " T=" + str(T))
		os.system("./glatzprogram triJumps=0 outpre=../../data/noTri/ writelines=500 Nt=" + str(Nt) + " Ex=" + str(Ex) + " Hz=0" + " T=" + str(T))
