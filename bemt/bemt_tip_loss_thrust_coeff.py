### Minor Project Code - 3
## Coded by : G R Krishna Chand Avatar, Kumar Gaurav, Kashish Korotania
## BEMT - Thrust coefficient Distribution and including the Tip Losses

import math
import numpy as np
import matplotlib.pyplot as plt
def prandtl_tip_loss(r,theta):
   global N_b, sol, cl_a, lam_c
   err = 1
   lam_old = (sol*cl_a/16)*(math.sqrt(1 + 32*theta*r/(sol*cl_a)) - 1)
   while(err > 1e-5):
        f = 0.5*N_b*(1-r)/lam_old
        F = (2/math.pi)*math.acos(math.exp(-f))
        lam_new = (sol*cl_a/(16*F))*(math.sqrt(1 + 32*F*theta*r/(sol*cl_a)) - 1)
        err = abs(lam_new - lam_old)
        lam_old = lam_new
   return(lam_old)
# Constant Parameters
sol = 0.0578 #Solidity
cl_a = 6.28 #Cl_a
N_b = 4   #Number of Blades
lam_c = 0
#Initialising Parameters                                                           
r = np.linspace(0,0.999999,500)
N = len(r)
# Computation begins
theta_o = 10   # pitch angle at root
theta_tw = -2.5 # linear twist rate
th = []     # Theta - pitch angle distribution
del_ct = [] # Coefficient of thrust
lam = [] # Inflow distribution
   
for i in range(N):
   theta = theta_o + theta_tw*r[i]
   th = th + [theta*math.pi/180]
for i in range(N):
   lam_t = prandtl_tip_loss(r[i],th[i])
   lam = lam + [lam_t]
   del_ct_t = N_b*0.5*(sol*cl_a)*(th[i]*r[i]**2 - lam_t*r[i])
   del_ct = del_ct + [del_ct_t]
# Writing to output file
f = "bemt_output_" + str(theta) + ".txt"
out = open(f,'w')
out.write('#r \t Inflow distribution\t\tThrust coefficient\n')
for i in range(N):
   out.write(str(round(r[i],9))+ '\t' + str(round(lam[i],9))+ '\t' + str(round(del_ct[i],9))+ '\n')
out.close()
# Plotting
plt.plot(r,del_ct,'b')
plt.xlabel('r')
plt.ylabel('dC_T/dr')
plt.xticks(np.arange(0,1,0.1))
plt.grid(True)
plt.title('Thrust coefficient distribution for theta_root = '+str(theta_o)+' degrees, linear twist rate = '+ str(theta_tw))
plt.show()
