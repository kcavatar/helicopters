## Minor Project code - 2
## Coded by : G R Krishna Chand Avatar, Kumar Gaurav, Kashish Korotania
## BEMT - Inflow Distribution and including the Tip Losses
import math
import numpy as np
import matplotlib.pyplot as plt
def prandtl_tip_loss(r,theta):
   global N_b, sol, cl_a, lam_c
   err = 1
   lam_old = (sol*cl_a/16)*(math.sqrt(1 + 32*theta*r/(sol*cl_a)) - 1)
   while(err > 1e-5):
        f = 0.5*N_b*(1-r)/lam_old 
		# f_root = 0.5*N_b*r**2/(1-r)/lam_old
        F = (2/math.pi)*math.acos(math.exp(-f))
        # lam_new = (sol*cl_a/(16*F))*(math.sqrt(1 + 32*F*theta*r/(sol*cl_a)) - 1)
        lam_new =  - (sol*cl_a/(16*F) - lam_c/2) + math.sqrt((sol*cl_a/16/F - lam_c/2)**2 + sol*cl_a*theta*r/8/F)
        err = abs(lam_new - lam_old)
        lam_old = lam_new
   return(lam_old)  
# Constant Parameters
sol = 0.0578 #Solidity
cl_a = 6.28 #Cl_a
N_b = 4   #Number of Blades
lam_c = 0.1  # Non-dimensional axial velocity
#Initialising Parameters                                                           
r = np.linspace(0,0.999999,500)
N = len(r)
# Computation for varying theta
theta_o = 10   # pitch angle at root
theta_tw = -2.5 # linear twist rate
th = []     # Theta - pitch angle distribution
del_ct = [] # Coefficient of thrust
lam = [] # Inflow distribution
lam_no_tip = []  
del_ct_no_tip = [] 
for i in range(N):
   theta = theta_o + theta_tw*r[i]
   th = th + [theta*math.pi/180]
for i in range(N):
   lam_t = prandtl_tip_loss(r[i],th[i])
   lam = lam + [lam_t]
   lam_no_tip_t = -(sol*cl_a/16 - lam_c/2) + math.sqrt((sol*cl_a/16 - lam_c/2)**2 + sol*cl_a*th[i]*r[i]/8)
   lam_no_tip = lam_no_tip + [lam_no_tip_t]
   del_ct_t = 0.5*(sol*cl_a)*(th[i]*r[i]**2 - lam_t*r[i])
   del_ct = del_ct + [del_ct_t]
   del_ct_no_tip_t = 0.5*(sol*cl_a)*(th[i]*r[i]**2 - lam_no_tip_t*r[i])
   del_ct_no_tip = del_ct_no_tip + [del_ct_no_tip_t]

# Plotting 
plt.plot(r,lam,'b', label='BEMT+Vortex theory')
plt.plot(r,lam_no_tip, 'm--', label='BEMT')
plt.xlabel('Non-dimensional radius (r)')
plt.ylabel('Inflow ratio (lambda)')
plt.legend(loc=4)
plt.xticks(np.arange(0,1.1,0.1))
#plt.ylabel('Inflow ratio (lam)')
plt.grid(True)
plt.title('Inflow ratio distribution for theta_root = '+str(theta_o)+' degrees, linear twist rate = '+ str(theta_tw))
plt.show()
