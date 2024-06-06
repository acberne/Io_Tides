
import numpy as np
import Get_Stress_and_Strain
import subprocess


#Prior to running this code, enter first enter the following command into shell:

#  ./tides.a elastic_model_Io.dat  #

#Then run this code. Specify co-latitude (thet), longitude (phi), and mean anomaly (M), all in degrees

thet= 0
phi=0
M=0

tau_thetthet, tau_phiphi, tau_phithet, eps_thetthet, eps_phiphi, eps_phithet = Get_Stress_and_Strain.Get_Stress(thet,phi,M)

np.savetxt('tau_thetthet.txt',[tau_thetthet])
np.savetxt('tau_phiphi.txt',[tau_phiphi])
np.savetxt('tau_phithet.txt',[tau_phithet])

np.savetxt('eps_thetthet.txt',[eps_thetthet])
np.savetxt('eps_phiphi.txt',[eps_phiphi])
np.savetxt('eps_phithet.txt',[eps_phithet])
